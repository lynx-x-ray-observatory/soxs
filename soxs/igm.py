from astropy.io import fits
from soxs.constants import K_per_keV
from soxs.utils import parse_value, get_data_file, mylog
import numpy as np


class IGMGenerator:
    def __init__(self, emin, emax, resonant_scattering=False, cxb_factor=1.0,
                 use_var_elem=False):
        self.binscale = "log"
        self.resonant_scattering = resonant_scattering
        self.cosmic_table = get_data_file("igm_v2ph_nome.fits")
        self.metal_tab_names = ["ox", "ne", "si", "su", "fe"] 
        if resonant_scattering:
            if use_var_elem:
                self.metal_tables = (get_data_file("igm_v2ph_mxx.fits"), 
                                     get_data_file("igm_v2sc_mxx.fits"))
                self.var_tables = [(get_data_file(f"igm_v2ph_{el}.fits"),
                                    get_data_file(f"igm_v2sc_{el}.fits"))
                                    for el in self.metal_tab_names]
            else:
                self.metal_tables = (get_data_file("igm_v2ph_me.fits"), 
                                     get_data_file("igm_v2sc_me.fits"))
        else:
            if use_var_elem:
                self.metal_tables = (get_data_file("igm_v2ph_mxx.fits"),)
                self.var_tables = [(get_data_file(f"igm_v2ph_{el}.fits"),)
                                   for el in self.metal_tab_names]
            else:
                self.metal_tables = (get_data_file("igm_v2ph_me.fits"),)
        self.cxb_factor = cxb_factor
        self.max_tables = 2 if resonant_scattering else 1
        self.var_elem = ["O", "Ne", "Si", "S", "Fe"] if use_var_elem else None
        self.nvar_elem = len(self.var_elem)
        self.emin = emin
        self.emax = emax
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
                for j, vfile in enumerate(self.var_tables):
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
            {"O": 0.4, "N": 0.3, "He": 0.9}. Default: None
        """
        from soxs.spectra import Spectrum
        if elem_abund is None:
            elem_abund = {}
        if set(elem_abund.keys()) != set(self.var_elem):
            raise RuntimeError("The supplied set of abundances is not the "
                               "same as that was originally set!\n"
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
                spec += dx1*vspec[j,idx1,:]
                spec += dx2*vspec[j,idx2,:]
                spec += dx3*vspec[j,idx3,:]
                spec += dx4*vspec[j,idx4,:]
        spec = 1.0e14*norm*spec[0,:]/de
        return Spectrum(ebins, spec, binscale=self.binscale)
