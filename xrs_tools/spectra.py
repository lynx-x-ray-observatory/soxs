import numpy as np
import subprocess
import tempfile
import shutil
import os

class Spectrum(object):
    def __init__(self, ebins, flux):
        self.ebins = ebins
        self.emid = 0.5*(ebins[1:]+ebins[:-1])
        self.flux = flux
        self.nbins = len(self.emid)
        self.tot_flux = self.flux.sum()

    @classmethod
    def from_xspec(cls, model_string, params, emin=0.01, emax=50.0,
                   nbins=10000):
        tmpdir = tempfile.mkdtemp()
        curdir = os.getcwd()
        os.chdir(tmpdir)
        xspec_in = []
        model_str = "%s &" % model_string
        for param in params:
            model_str += " %g &" % param
        model_str += " /*"
        xspec_in.append("model %s\n" % model_str)
        xspec_in.append("dummyrsp %g %g %d lin\n" % (emin, emax, nbins))
        xspec_in += ["set fp [open spec_therm.xspec w+]\n",
                     "tclout energies\n", "puts $fp $xspec_tclout\n",
                     "tclout modval\n", "puts $fp $xspec_tclout\n",
                     "close $fp\n", "quit\n"]
        f_xin = open("xspec.in", "w")
        f_xin.writelines(xspec_in)
        f_xin.close()
        subprocess.call(["xspec", "-", "xspec.in"])
        f_s = open("spec_therm.xspec", "r")
        lines = f_s.readlines()
        f_s.close()
        ebins = np.array(lines[0].split()).astype("float64")
        flux = np.array(lines[1].split()).astype("float64")
        os.chdir(curdir)
        shutil.rmtree(tmpdir)
        return cls(ebins, flux)

    @classmethod
    def from_apec(cls, kT, abund, redshift, norm, 
                  broadening=False, absorb_model='tbabs',
                  nH=0.02, velocity=0.0, emin=0.01, 
                  emax=50.0, nbins=10000):
        if broadening:
            model_str = "bapec"
            params = [kT, abund, redshift, velocity, norm]
        else:
            model_str = "apec"
            params = [kT, abund, redshift, norm]
        if absorb_model is not None:
            model_str = '*'.join([absorb_model, model_str])
            params = [nH] + params
        return cls.from_xspec(model_str, params, emin=emin, emax=emax,
                              nbins=nbins)

    @classmethod
    def from_powerlaw(cls, photon_index, redshift, norm,
                      absorb_model='tbabs', nH=0.02,
                      emin=0.01, emax=50.0, nbins=10000):
        model_str = "zpowerlw"
        params = [photon_index, redshift, norm]
        if absorb_model is not None:
            model_str = '*'.join([absorb_model, model_str])
            params = [nH] + params
        return cls.from_xspec(model_str, params, emin=emin, emax=emax,
                              nbins=nbins)

    @classmethod
    def from_file(cls, filename):
        data = np.loadtxt(filename)
        emid = data[:,0]
        if data.shape[-1] == 3:
            de = data[:,1]
        else:
            de = np.diff(emid)[0]
        flux = data[:,-1]
        ebins = np.concatenate([emid-0.5*de, emid[-1]+0.5*de])
        return cls(ebins, flux)

    def rescale_flux(self, new_flux, emin=None, emax=None):
        """
        Rescale the flux of the spectrum, optionally using a 
        specific energy band.

        Parameters
        ----------
        new_flux : float
            The new flux in units of photons/s/cm**2.
        emin : float, optional
            The minimum energy of the band to consider, in keV. 
            Default: Use the minimum energy of the entire spectrum.
        emax : float, optional
            The maximum energy of the band to consider, in keV.
            Default: Use the maximum energy of the entire spectrum.
        """
        if emin is None:
            emin = self.ebins[0]
        if emax is None:
            emax = self.ebins[-1]
        idxs = np.logical_and(self.emid >= emin, self.emid <= emax)
        f = self.flux[idxs].sum()
        self.flux *= new_flux/f
        self.tot_flux = self.flux.sum()

    def generate_energies(self, t_exp, area, prng=None):
        """
        Generate photon energies from this spectrum given an 
        exposure time and effective area. 

        Parameters
        ----------
        t_exp : float
            The exposure time in seconds. 
        area : float
            The effective area (constant) in cm**2.
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        """
        if prng is None:
            prng = np.random
        cumspec = np.cumsum(self.flux)
        n_ph = t_exp*area*cumspec[-1]
        n_ph = np.uint64(n_ph) + np.uint64(n_ph >= prng.uniform())
        cumspec = np.insert(cumspec, 0, 0.0)
        cumspec /= cumspec[-1]
        randvec = prng.uniform(size=n_ph)
        randvec.sort()
        energies = np.interp(randvec, cumspec, self.ebins)
        return energies

def determine_apec_norm(kT, abund, redshift, emin, emax,
                        flux, nbins=10000):
    """
    Determine the normalization of an unabsorbed APEC 
    model given a total flux within some band.

    Parameters
    ----------
    kT : float
        The temperature of the model in keV.
    abund : float
        The abundance in solar units.
    redshift : float
        The redshift.
    emin : float
        The lower energy of the energy band to consider, in keV.
    emax : float
        The upper energy of the energy band to consider, in keV.
    flux : float
        The total flux in the band in photons/s/cm**2. 
    nbins : integer, optional
        The number of bins to sum the spectrum over. 
        Default: 10000
    """
    spec = Spectrum.from_apec(kT, abund, redshift, 1.0, 
                              emin=emin, emax=emax, nbins=nbins,
                              absorb_model=None)
    return flux/spec.tot_flux

def determine_powerlaw_norm(photon_index, redshift, emin, emax,
                            flux, nbins=10000):
    """
    Determine the normalization of an unabsorbed APEC 
    model given a total flux within some band.

    Parameters
    ----------
    photon_index : float
        The temperature of the model in keV.
    redshift : float
        The redshift.
    emin : float
        The lower energy of the energy band to consider, in keV.
    emax : float
        The upper energy of the energy band to consider, in keV.
    flux : float
        The total flux in the band in photons/s/cm**2. 
    nbins : integer, optional
        The number of bins to sum the spectrum over. 
        Default: 10000
    """
    spec = Spectrum.from_apec(photon_index, redshift, 1.0,
                              emin=emin, emax=emax, nbins=nbins,
                              absorb_model=None)
    return flux/spec.tot_flux
