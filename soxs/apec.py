import numpy as np
from astropy.io import fits
import os
from soxs.utils import parse_value, mylog, soxs_cfg, \
    DummyPbar, get_data_file
from tqdm.auto import tqdm
from soxs.lib.broaden_lines import broaden_lines
from soxs.constants import erg_per_keV, hc, \
    cosmic_elem, metal_elem, atomic_weights, clight, \
    m_u, elem_names, abund_tables


class ApecGenerator:
    r"""
    Initialize a thermal gas emission model from the 
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
        "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914 
        except for elements not listed which are given zero abundance)
        "lodd" : from Lodders, K (2003, ApJ 591, 1220)
    nei : boolean, optional
        If True, use the non-equilibrium ionization tables. These are
        not supplied with SOXS but must be downloaded separately, in
        which case the *apec_root* parameter must also be set to their
        location. Default: False

    Examples
    --------
    >>> apec_model = ApecGenerator(0.05, 50.0, 1000, apec_vers="3.0.3",
    ...                            broadening=True)
    """
    def __init__(self, emin, emax, nbins, var_elem=None, apec_root=None,
                 apec_vers=None, broadening=True, nolines=False,
                 abund_table=None, nei=False):
        if apec_vers is None:
            apec_vers = soxs_cfg.get("soxs", "apec_vers")
        mylog.debug(f"Using APEC version {apec_vers}.")
        if nei and var_elem is None:
            raise RuntimeError("For NEI spectra, you must specify which elements "
                               "you want to vary using the 'var_elem' argument!")
        self.nei = nei
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, 'keV')
        self.emin = emin
        self.emax = emax
        self.nbins = nbins
        self.ebins = np.linspace(self.emin, self.emax, nbins+1)
        self.de = np.diff(self.ebins)
        self.emid = 0.5*(self.ebins[1:]+self.ebins[:-1])
        if nei:
            neistr = "_nei"
            ftype = "comp"
        else:
            neistr = ""
            ftype = "coco"
        cocofile = f"apec_v{apec_vers}{neistr}_{ftype}.fits"
        linefile = f"apec_v{apec_vers}{neistr}_line.fits"
        if apec_root is None:
            self.cocofile = get_data_file(cocofile)
            self.linefile = get_data_file(linefile)
        else:
            self.cocofile = os.path.join(apec_root, cocofile)
            self.linefile = os.path.join(apec_root, linefile)
        if not os.path.exists(self.cocofile) or not os.path.exists(self.linefile):
            raise IOError(f"Cannot find the APEC files!\n {self.cocofile}\n, "
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

    def _spectrum_init(self, kT, velocity, elem_abund):
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
        kT, dT, tindex, v = self._spectrum_init(kT, velocity, elem_abund)
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return np.zeros(self.nbins)
        cspec, mspec, vspec = self._get_table([tindex, tindex+1], redshift, v)
        cosmic_spec = cspec[0,:]*(1.-dT)+cspec[1,:]*dT
        metal_spec = mspec[0,:]*(1.-dT)+mspec[1,:]*dT
        spec = cosmic_spec + abund*metal_spec
        if vspec is not None:
            for elem, eabund in elem_abund.items():
                j = self.var_elem_names.index(elem)
                spec += eabund*(vspec[j,0,:]*(1.-dT)+vspec[j,1,:]*dT)
        spec = 1.0e14*norm*spec/self.de
        return Spectrum(self.ebins, spec)

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
        kT, dT, tindex, v = self._spectrum_init(kT, velocity, elem_abund)
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return np.zeros(self.nbins)
        cspec, _, vspec = self._get_table([tindex, tindex+1], redshift, v)
        spec = cspec[0,:]*(1.-dT)+cspec[1,:]*dT
        for elem, eabund in elem_abund.items():
            j = self.var_ion_names.index(elem)
            spec += eabund*(vspec[j,0,:]*(1.-dT) + vspec[j,1,:]*dT)
        spec = 1.0e14*norm*spec/self.de
        return Spectrum(self.ebins, spec)
