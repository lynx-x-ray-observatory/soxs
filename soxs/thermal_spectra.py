import numpy as np
from astropy.io import fits
import os
from soxs.utils import parse_value, mylog, soxs_cfg, \
    DummyPbar, get_data_file
from tqdm.auto import tqdm
from soxs.lib.broaden_lines import broaden_lines
from soxs.constants import erg_per_keV, hc, \
    cosmic_elem, metal_elem, atomic_weights, clight, \
    m_u, elem_names, abund_tables, K_per_keV


class CIEGenerator:
    def __init__(self, model, emin, emax, nbins, binscale='linear', 
                 var_elem=None, model_root=None,
                 model_vers=None, broadening=True, nolines=False,
                 abund_table=None, nei=False):
        self.binscale = binscale
        if model_vers is None:
            model_vers = soxs_cfg.get("soxs", f"{model.lower()}_vers")
        mylog.debug(f"Using {model} version {model_vers}.")
        if nei and var_elem is None:
            raise RuntimeError("For NEI spectra, you must specify which elements "
                               "you want to vary using the 'var_elem' argument!")
        self.nei = nei
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, 'keV')
        self.emin = emin
        self.emax = emax
        self.nbins = nbins
        if binscale == "linear":
            self.ebins = np.linspace(self.emin, self.emax, nbins+1)
        elif binscale == "log":
            self.ebins = np.logspace(np.log10(self.emin), np.log10(self.emax), nbins+1)
        self.de = np.diff(self.ebins)
        self.emid = 0.5*(self.ebins[1:]+self.ebins[:-1])
        if nei:
            neistr = "_nei"
            ftype = "comp"
        else:
            neistr = ""
            ftype = "coco"
        cocofile = f"{model.lower()}_v{model_vers}{neistr}_{ftype}.fits"
        linefile = f"{model.lower()}_v{model_vers}{neistr}_line.fits"
        if model_root is None:
            self.cocofile = get_data_file(cocofile)
            self.linefile = get_data_file(linefile)
        else:
            self.cocofile = os.path.join(model_root, cocofile)
            self.linefile = os.path.join(model_root, linefile)
        if not os.path.exists(self.cocofile) or not os.path.exists(self.linefile):
            raise IOError(f"Cannot find the {model} files!\n {self.cocofile}\n, "
                          f"{self.linefile}")
        mylog.debug(f"Using {cocofile} for generating spectral lines.")
        mylog.debug(f"Using {linefile} for generating the continuum.")
        self.nolines = nolines
        self.wvbins = hc/self.ebins[::-1]
        self.broadening = broadening
        self.line_handle = fits.open(self.linefile)
        self.coco_handle = fits.open(self.cocofile)
        self.nT = self.line_handle[1].data.shape[0]
        self.Tvals = self.line_handle[1].data.field("kT")
        self.dTvals = np.diff(self.Tvals)
        self.minlam = self.wvbins.min()
        self.maxlam = self.wvbins.max()
        self.var_elem_names = []
        self.var_ion_names = []
        if var_elem is None:
            self.var_elem = np.empty((0, 1), dtype='int')
        else:
            self.var_elem = []
            if len(var_elem) != len(set(var_elem)):
                raise RuntimeError("Duplicates were found in the \"var_elem\" list! %s" % var_elem)
            for elem in var_elem:
                if "^" in elem:
                    if not self.nei:
                        raise RuntimeError("Cannot use different ionization states with a "
                                           "CIE plasma!")
                    el = elem.split("^")
                    e = el[0]
                    ion = int(el[1])
                else:
                    if self.nei:
                        raise RuntimeError("Variable elements must include the ionization "
                                           "state for NEI plasmas!")
                    e = elem
                    ion = 0
                self.var_elem.append([elem_names.index(e), ion])
            self.var_elem.sort(key=lambda x: (x[0], x[1]))
            self.var_elem = np.array(self.var_elem, dtype='int')
            self.var_elem_names = [elem_names[e[0]] for e in self.var_elem]
            self.var_ion_names = ["%s^%d" % (elem_names[e[0]], e[1]) for e in self.var_elem]
        self.num_var_elem = len(self.var_elem)
        if self.nei:
            self.cosmic_elem = [elem for elem in [1, 2]
                                if elem not in self.var_elem[:, 0]]
            self.metal_elem = []
        else:
            self.cosmic_elem = [elem for elem in cosmic_elem
                                if elem not in self.var_elem[:,0]]
            self.metal_elem = [elem for elem in metal_elem
                               if elem not in self.var_elem[:,0]]
        if abund_table is None:
            abund_table = soxs_cfg.get("soxs", "abund_table")
        if not isinstance(abund_table, str):
            if len(abund_table) != 30:
                raise RuntimeError("User-supplied abundance tables "
                                   "must be 30 elements long!")
            self.atable = np.concatenate([[0.0], np.array(abund_table)])
        else:
            self.atable = abund_tables[abund_table].copy()
        self._atable = self.atable.copy()
        self._atable[1:] /= abund_tables["angr"][1:]

    def _make_spectrum(self, kT, element, ion, velocity, line_fields,
                       coco_fields, scale_factor):

        tmpspec = np.zeros(self.nbins)

        if not self.nolines:
            loc = (line_fields['element'] == element) & \
                  (line_fields['lambda'] > self.minlam) & \
                  (line_fields['lambda'] < self.maxlam)
            if self.nei:
                loc &= (line_fields['ion_drv'] == ion+1)
            i = np.where(loc)[0]
            E0 = hc/line_fields['lambda'][i].astype("float64")*scale_factor
            amp = line_fields['epsilon'][i].astype("float64")*self._atable[element]
            if self.broadening:
                sigma = 2.*kT*erg_per_keV/(atomic_weights[element]*m_u)
                sigma += 2.0*velocity*velocity
                sigma = E0*np.sqrt(sigma)/clight
                vec = broaden_lines(E0, sigma, amp, self.ebins)
            else:
                vec = np.histogram(E0, self.ebins, weights=amp)[0]
            tmpspec += vec

        ind = np.where((coco_fields['Z'] == element) &
                       (coco_fields['rmJ'] == ion+int(self.nei)))[0]

        if len(ind) == 0:
            return tmpspec
        else:
            ind = ind[0]

        de0 = self.de/scale_factor

        n_cont = coco_fields['N_Cont'][ind]
        e_cont = coco_fields['E_Cont'][ind][:n_cont]*scale_factor
        continuum = coco_fields['Continuum'][ind][:n_cont]*self._atable[element]

        tmpspec += np.interp(self.emid, e_cont, continuum)*de0

        n_pseudo = coco_fields['N_Pseudo'][ind]
        e_pseudo = coco_fields['E_Pseudo'][ind][:n_pseudo]*scale_factor
        pseudo = coco_fields['Pseudo'][ind][:n_pseudo]*self._atable[element]

        tmpspec += np.interp(self.emid, e_pseudo, pseudo)*de0

        return tmpspec*scale_factor

    def _preload_data(self, index):
        line_data = self.line_handle[index+2].data
        coco_data = self.coco_handle[index+2].data
        line_fields = ['element', 'lambda', 'epsilon']
        if self.nei:
            line_fields.append('ion_drv')
        line_fields = tuple(line_fields)
        coco_fields = ('Z', 'rmJ', 'N_Cont', 'E_Cont', 'Continuum',
                       'N_Pseudo','E_Pseudo', 'Pseudo')
        line_fields = {el: line_data.field(el) for el in line_fields}
        coco_fields = {el: coco_data.field(el) for el in coco_fields}
        return line_fields, coco_fields

    def _get_table(self, indices, redshift, velocity):
        numi = len(indices)
        scale_factor = 1./(1.+redshift)
        cspec = np.zeros((numi, self.nbins))
        mspec = np.zeros((numi, self.nbins))
        vspec = None
        if self.num_var_elem > 0:
            vspec = np.zeros((self.num_var_elem, numi, self.nbins))
        if numi > 2:
            pbar = tqdm(leave=True, total=numi, desc="Preparing spectrum table ")
        else:
            pbar = DummyPbar()
        for i, ikT in enumerate(indices):
            line_fields, coco_fields = self._preload_data(ikT)
            # First do H, He, and trace elements
            for elem in self.cosmic_elem:
                if self.nei:
                    # For H, He we assume fully ionized
                    ion = elem
                else:
                    ion = 0
                cspec[i,:] += self._make_spectrum(self.Tvals[ikT], elem, ion, velocity, line_fields,
                                                  coco_fields, scale_factor)
            # Next do the metals
            for elem in self.metal_elem:
                mspec[i,:] += self._make_spectrum(self.Tvals[ikT], elem, 0, velocity, line_fields,
                                                  coco_fields, scale_factor)
            # Now do any metals that we wanted to vary freely from the abund
            # parameter
            if self.num_var_elem > 0:
                for j, elem in enumerate(self.var_elem):
                    vspec[j,i,:] = self._make_spectrum(self.Tvals[ikT], elem[0], elem[1],
                                                       velocity, line_fields, coco_fields, scale_factor)
            pbar.update()
        pbar.close()
        return cspec, mspec, vspec

    def _spectrum_init(self, kT, velocity):
        kT = parse_value(kT, "keV")
        velocity = parse_value(velocity, "km/s")
        v = velocity*1.0e5
        tindex = np.searchsorted(self.Tvals, kT)-1
        dT = (kT-self.Tvals[tindex])/self.dTvals[tindex]
        return kT, dT, tindex, v

    def get_spectrum(self, kT, abund, redshift, norm, velocity=0.0,
                     elem_abund=None):
        """
        Get a thermal emission spectrum assuming CIE.

        Parameters
        ----------
        kT : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The temperature in keV.
        abund : float
            The metal abundance in solar units. 
        redshift : float
            The redshift.
        norm : float
            The normalization of the model, in the standard
            Xspec units of 1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
        velocity : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The velocity broadening parameter, in units of 
            km/s. Default: 0.0
        elem_abund : dict of element name, float pairs, optional
            A dictionary of elemental abundances in solar
            units to vary freely of the abund parameter, e.g.
            {"O": 0.4, "N": 0.3, "He": 0.9}. Default: None
        """
        from soxs.spectra import Spectrum
        if self.nei:
            raise RuntimeError("Use 'get_nei_spectrum' for NEI spectra!")
        if elem_abund is None:
            elem_abund = {}
        if set(elem_abund.keys()) != set(self.var_elem_names):
            raise RuntimeError("The supplied set of abundances is not the "
                               "same as that was originally set!\n"
                               "Free elements: %s\nAbundances: %s" % (set(elem_abund.keys()),
                                                                      set(self.var_elem_names)))
        kT, dT, tindex, v = self._spectrum_init(kT, velocity)
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return Spectrum(self.ebins, np.zeros(self.nbins), binscale=self.binscale)
        cspec, mspec, vspec = self._get_table([tindex, tindex+1], redshift, v)
        cosmic_spec = cspec[0,:]*(1.-dT)+cspec[1,:]*dT
        metal_spec = mspec[0,:]*(1.-dT)+mspec[1,:]*dT
        spec = cosmic_spec + abund*metal_spec
        if vspec is not None:
            for elem, eabund in elem_abund.items():
                j = self.var_elem_names.index(elem)
                spec += eabund*(vspec[j,0,:]*(1.-dT)+vspec[j,1,:]*dT)
        spec = 1.0e14*norm*spec/self.de
        return Spectrum(self.ebins, spec, binscale=self.binscale)

    def get_nei_spectrum(self, kT, elem_abund, redshift, norm, velocity=0.0):
        """
        Get a thermal emission spectrum assuming NEI.

        Parameters
        ----------
        kT : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The temperature in keV.
        elem_abund : dict of element name, float pairs
            A dictionary of ionization state abundances in solar
            units to vary freely of the abund parameter, e.g.
            {"O^1": 0.4, "O^4": 0.6, "N^2": 0.7} Default: None
        redshift : float
            The redshift.
        norm : float
            The normalization of the model, in the standard
            Xspec units of 1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
        velocity : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The velocity broadening parameter, in units of 
            km/s. Default: 0.0
        """
        from soxs.spectra import Spectrum
        if not self.nei:
            raise RuntimeError("Use 'get_spectrum' for CIE spectra!")
        if set(elem_abund.keys()) != set(self.var_ion_names):
            raise RuntimeError("The supplied set of abundances is not the "
                               "same as that was originally set!\n"
                               "Free elements: %s\nAbundances: %s" % (set(elem_abund.keys()),
                                                                      set(self.var_ion_names)))
        kT, dT, tindex, v = self._spectrum_init(kT, velocity)
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return np.zeros(self.nbins)
        cspec, _, vspec = self._get_table([tindex, tindex+1], redshift, v)
        spec = cspec[0,:]*(1.-dT)+cspec[1,:]*dT
        for elem, eabund in elem_abund.items():
            j = self.var_ion_names.index(elem)
            spec += eabund*(vspec[j,0,:]*(1.-dT) + vspec[j,1,:]*dT)
        spec = 1.0e14*norm*spec/self.de
        return Spectrum(self.ebins, spec, binscale=self.binscale)


class ApecGenerator(CIEGenerator):
    r"""
    Create spectra assuming a thermal plasma emission model
    in collisional ionization equilibrium from the 
    AtomDB APEC tables available at http://www.atomdb.org. 
    This code borrows heavily from Python routines used to 
    read the APEC tables developed by Adam Foster at the
    CfA (afoster@cfa.harvard.edu).

    Parameters
    ----------
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The minimum energy for the spectral model.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The maximum energy for the spectral model.
    nbins : integer
        The number of bins in the spectral model.
    binscale : string, optional
        The scale of the energy binning: "linear" or "log". 
        Default: "linear"
    var_elem : list of strings, optional
        The names of elements to allow to vary freely
        from the single abundance parameter. These can be strings like
        ["O", "N", "He"], or if nei=True they must be elements with
        ionization states, e.g. ["O^1", "O^2", "N^4"]. Default:
        None
    apec_root : string, optional
        The directory root where the APEC model files 
        are stored. If not provided, the default is to 
        grab them from the tables stored with SOXS.
    apec_vers : string, optional
        The version identifier string for the APEC files. Default is
        set in the SOXS configuration file, the default for which is
        "3.0.9".
    broadening : boolean, optional
        Whether or not the spectral lines should be 
        thermally and velocity broadened. Default: True
    nolines : boolean, optional
        Turn off lines entirely for generating spectra.
        Default: False
    abund_table : string or array_like, optional
        The abundance table to be used for solar abundances. 
        Either a string corresponding to a built-in table or an array
        of 30 floats corresponding to the abundances of each element
        relative to the abundance of H. Default is set in the SOXS
        configuration file, the default for which is "angr".
        Built-in options are:
        "angr" : from Anders E. & Grevesse N. (1989, Geochimica et 
        Cosmochimica Acta 53, 197)
        "aspl" : from Asplund M., Grevesse N., Sauval A.J. & Scott 
        P. (2009, ARAA, 47, 481)
        "feld" : from Feldman U. (1992, Physica Scripta, 46, 202)
        "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914 
        except for elements not listed which are given zero abundance)
        "lodd" : from Lodders, K (2003, ApJ 591, 1220)
    nei : boolean, optional
        If True, use the non-equilibrium ionization tables.

    Examples
    --------
    >>> apec_model = ApecGenerator(0.05, 50.0, 1000, apec_vers="3.0.3",
    ...                            broadening=True)
    """
    def __init__(self, emin, emax, nbins, binscale='linear',
                 var_elem=None, apec_root=None,
                 apec_vers=None, broadening=True, nolines=False,
                 abund_table=None, nei=False):
        super().__init__("apec", emin, emax, nbins, binscale=binscale, 
                         var_elem=var_elem, model_root=apec_root,
                         model_vers=apec_vers, broadening=broadening,
                         nolines=nolines, abund_table=abund_table,
                         nei=nei)


class SpexGenerator(CIEGenerator):
    r"""
    Create thermal APEC tables available at http://www.atomdb.org. 
    The same underlying machinery as the APEC model is used, as the
    SPEX model has been converted to the APEC table format using the
    code at https://github.com/jeremysanders/spex_to_xspec.

    Parameters
    ----------
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The minimum energy for the spectral model.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The maximum energy for the spectral model.
    nbins : integer
        The number of bins in the spectral model.
    binscale : string, optional
        The scale of the energy binning: "linear" or "log". 
        Default: "linear"
    var_elem : list of strings, optional
        The names of elements to allow to vary freely
        from the single abundance parameter. These can be strings like
        ["O", "N", "He"], or if nei=True they must be elements with
        ionization states, e.g. ["O^1", "O^2", "N^4"]. Default:
        None
    spex_root : string, optional
        The directory root where the SPEX model files 
        are stored. If not provided, the default is to 
        grab them from the tables stored with SOXS.
    spex_vers : string, optional
        The version identifier string for the SPEX files. Default is
        set in the SOXS configuration file, the default for which is
        "3.06.01".
    broadening : boolean, optional
        Whether or not the spectral lines should be 
        thermally and velocity broadened. Default: True
    nolines : boolean, optional
        Turn off lines entirely for generating spectra.
        Default: False
    abund_table : string or array_like, optional
        The abundance table to be used for solar abundances. 
        Either a string corresponding to a built-in table or an array
        of 30 floats corresponding to the abundances of each element
        relative to the abundance of H. Default is set in the SOXS
        configuration file, the default for which is "angr".
        Built-in options are:
        "angr" : from Anders E. & Grevesse N. (1989, Geochimica et 
        Cosmochimica Acta 53, 197)
        "aspl" : from Asplund M., Grevesse N., Sauval A.J. & Scott 
        P. (2009, ARAA, 47, 481)
        "feld" : from Feldman U. (1992, Physica Scripta, 46, 202)
        "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914 
        except for elements not listed which are given zero abundance)
        "lodd" : from Lodders, K (2003, ApJ 591, 1220)

    Examples
    --------
    >>> spex_model = SpexGenerator(0.05, 50.0, 1000, broadening=True)
    """
    def __init__(self, emin, emax, nbins, binscale='linear',
                 var_elem=None, spex_vers=None, spex_root=None,
                 broadening=True, nolines=False,
                 abund_table=None):
        super().__init__("spex", emin, emax, nbins, binscale=binscale,
                         var_elem=var_elem, model_root=spex_root,
                         model_vers=spex_vers, broadening=broadening,
                         nolines=nolines, abund_table=abund_table,
                         nei=False)

    def get_nei_spectrum(self, kT, elem_abund, redshift, norm, velocity=0.0):
        raise NotImplementedError


metal_tab_names = {
    "O": "ox",
    "Ne": "ne",
    "Mg": "mg",
    "Si": "si",
    "S": "su",
    "Fe": "fe"
}


class IGMGenerator:
    def __init__(self, emin, emax, resonant_scattering=False, cxb_factor=0.5,
                 var_elem_option=None):
        """
        Initialize an emission model for a thermal plasma including 
        photoionization and resonant scattering from the CXB based on 
        Khabibullin & Churazov 2019
        (https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4972K/) and Churazov 
        et al. 2001 (https://ui.adsabs.harvard.edu/abs/2001MNRAS.323...93C/).

        Assumes the abundance tables from Feldman 1992.

        Table data and README files can be found at
        https://wwwmpa.mpa-garching.mpg.de/~ildar/igm/v2x/.

        Parameters
        ----------
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy for the spectral model.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy for the spectral model.
        resonant_scattering : boolean, optional
            Whether or not to include the effects of resonant scattering
            from CXB photons. Default: False
        cxb_factor : float, optional
            The fraction of the CXB photons that are resonant scattered to enhance
            the lines. Default: 0.5
        var_elem_option: integer, optional
            An integer to choose between options for variable elements, which are:
            1: specify abundances of O, Ne, and Fe separately from other metals
            2: specify abundances of O, Ne, Mg, Si, S, and Fe separately from other
               metals
            Default: None, which means no metal abundances can be specified
            separately.
        """
        if var_elem_option is None:
            metal_option = "me"
            self.metal_tab_names = []
            self.var_elem = []
        elif var_elem_option == 1:
            metal_option = "mx"
            self.var_elem = ["O", "Ne", "Fe"]
        elif var_elem_option == 2:
            metal_option = "mxx"
            self.var_elem = ["O", "Ne", "Mg", "Si", "S", "Fe"]
        else:
            raise RuntimeError(f"Unsupported 'var_elem_option' = {var_elem_option}!")
        self.var_elem_option = var_elem_option
        self.binscale = "log"
        self.resonant_scattering = resonant_scattering
        self.cosmic_table = get_data_file("igm_v2ph2_nome.fits")
        if resonant_scattering:
            self.metal_tables = (get_data_file(f"igm_v2ph2_{metal_option}.fits"),
                                 get_data_file(f"igm_v2sc_{metal_option}.fits"))
            if var_elem_option:
                self.var_tables = [(get_data_file(f"igm_v2ph2_{metal_tab_names[el]}.fits"),
                                    get_data_file(f"igm_v2sc_{metal_tab_names[el]}.fits"))
                                   for el in self.var_elem]
        else:
            self.metal_tables = (get_data_file(f"igm_v2ph2_{metal_option}.fits"),)
            if var_elem_option:
                self.var_tables = [(get_data_file(f"igm_v2ph2_{metal_tab_names[el]}.fits"),)
                                   for el in self.var_elem]
        self.cxb_factor = cxb_factor
        self.nvar_elem = len(self.var_elem)
        self.emin = parse_value(emin, "keV")
        self.emax = parse_value(emax, "keV")
        with fits.open(self.cosmic_table) as f:
            self.n_D = f["PARAMETERS"].data["NUMBVALS"][0]
            self.Dvals = f["PARAMETERS"].data["VALUE"][0][:self.n_D]
            self.n_T = f["PARAMETERS"].data["NUMBVALS"][1]
            self.Tvals = f["PARAMETERS"].data["VALUE"][1][:self.n_T]
        self.dDvals = np.diff(self.Dvals)
        self.dTvals = np.diff(self.Tvals)

    def _get_energies(self, redshift):
        scale_factor = 1.0/(1.0+redshift)
        with fits.open(self.cosmic_table) as f:
            elo = f["ENERGIES"].data["ENERG_LO"]*scale_factor
            ehi = f["ENERGIES"].data["ENERG_HI"]*scale_factor
        eidxs = elo > self.emin
        eidxs &= ehi < self.emax
        ne = eidxs.sum()
        ebins = np.append(elo[eidxs], ehi[eidxs][-1])
        emid = 0.5*(ebins[1:]+ebins[:-1])
        de = np.diff(ebins)
        return eidxs, ne, ebins, emid, de

    def _get_table(self, ne, eidxs, redshift):
        norm_fac = 5.50964e-5*np.array([1.0, self.cxb_factor])
        scale_factor = 1.0/(1.0+redshift)
        metal_spec = np.zeros((self.n_T*self.n_D, ne))
        var_spec = None
        mylog.debug(f"Opening {self.cosmic_table}.")
        with fits.open(self.cosmic_table) as f:
            cosmic_spec = f["SPECTRA"].data["INTPSPEC"][:,eidxs]*norm_fac[0]
        cosmic_spec *= scale_factor
        for j, mfile in enumerate(self.metal_tables):
            mylog.debug(f"Opening {mfile}.")
            with fits.open(mfile) as f:
                metal_spec += f["SPECTRA"].data["INTPSPEC"][:,eidxs]*norm_fac[j]
        metal_spec *= scale_factor
        if self.nvar_elem > 0:
            var_spec = np.zeros((self.nvar_elem,self.n_T*self.n_D,ne))
            for i in range(self.nvar_elem):
                for j, vfile in enumerate(self.var_tables[i]):
                    mylog.debug(f"Opening {vfile}.")
                    with fits.open(vfile) as f:
                        var_spec[i,:,:] += f["SPECTRA"].data["INTPSPEC"][:, eidxs]*norm_fac[j]
            var_spec *= scale_factor
        return cosmic_spec, metal_spec, var_spec

    def get_spectrum(self, kT, nH, abund, redshift, norm, elem_abund=None):
        """
        Get an emission spectrum from the IGM model.

        Parameters
        ----------
        kT : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The temperature in keV.
        nH : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The hydrogen number density in cm**-3.
        abund : float
            The metal abundance in solar units. 
        redshift : float
            The redshift.
        norm : float
            The normalization of the model, in the standard
            Xspec units of 1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
        elem_abund : dict of element name, float pairs, optional
            A dictionary of elemental abundances in solar
            units to vary freely of the abund parameter, e.g.
            {"O": 0.4, "Ne": 0.3, "Fe": 0.9}. Default: None
        """
        from soxs.spectra import Spectrum
        if elem_abund is None:
            elem_abund = {}
        if set(elem_abund.keys()) != set(self.var_elem):
            raise RuntimeError("The supplied set of abundances does not match "
                               "what is available for 'var_elem_option = "
                               f"{self.var_elem_option}!\n"
                               "Free elements: %s\nAbundances: %s" % (set(elem_abund.keys()),
                                                                      set(self.var_elem)))
        kT = parse_value(kT, "keV")
        nH = parse_value(nH, "cm**-3")
        lkT = np.atleast_1d(np.log10(kT*K_per_keV))
        lnH = np.atleast_1d(np.log10(nH))
        tidx = np.searchsorted(self.Tvals, lkT)-1
        didx = np.searchsorted(self.Dvals, lnH)-1
        eidxs, ne, ebins, emid, de = self._get_energies(redshift)
        if tidx >= self.Tvals.size-1 or tidx < 0:
            return Spectrum(self.ebins, np.zeros(ne), binscale=self.binscale)
        if didx >= self.Dvals.size-1 or didx < 0:
            return Spectrum(self.ebins, np.zeros(ne), binscale=self.binscale)
        cspec, mspec, vspec = self._get_table(ne, eidxs, redshift)
        dT = (lkT - self.Tvals[tidx]) / self.dTvals[tidx]
        dn = (lnH - self.Dvals[didx]) / self.dDvals[didx]
        idx1 = np.ravel_multi_index((didx+1,tidx+1), (self.n_D, self.n_T))
        idx2 = np.ravel_multi_index((didx+1,tidx), (self.n_D, self.n_T))
        idx3 = np.ravel_multi_index((didx,tidx+1), (self.n_D, self.n_T))
        idx4 = np.ravel_multi_index((didx,tidx), (self.n_D, self.n_T))
        dx1 = dT*dn
        dx2 = dn-dx1
        dx3 = dT-dx1
        dx4 = 1.0+dx1-dT-dn
        cosmic_spec = dx1*cspec[idx1,:]
        cosmic_spec += dx2*cspec[idx2,:]
        cosmic_spec += dx3*cspec[idx3,:]
        cosmic_spec += dx4*cspec[idx4,:]
        metal_spec = dx1*mspec[idx1,:]
        metal_spec += dx2*mspec[idx2,:]
        metal_spec += dx3*mspec[idx3,:]
        metal_spec += dx4*mspec[idx4,:]
        spec = cosmic_spec + abund*metal_spec
        if vspec is not None:
            for elem, eabund in elem_abund.items():
                j = self.var_elem.index(elem)
                spec += eabund*dx1*vspec[j,idx1,:]
                spec += eabund*dx2*vspec[j,idx2,:]
                spec += eabund*dx3*vspec[j,idx3,:]
                spec += eabund*dx4*vspec[j,idx4,:]
        spec = norm*spec[0,:]/de/nH
        return Spectrum(ebins, spec, binscale=self.binscale)
