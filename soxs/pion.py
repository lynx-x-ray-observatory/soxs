from astropy.io import fits
from soxs.constants import K_per_keV
from soxs.utils import parse_value, get_data_file, mylog
import numpy as np


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
