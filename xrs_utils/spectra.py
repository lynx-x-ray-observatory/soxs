import numpy as np

class Spectrum(object):
    def __init__(self, ebins, flux):
        self.ebins = ebins
        self.emid = 0.5*(ebins[1:]+ebins[:-1])
        self.flux = flux
        self.nbins = len(self.emid)
        self.tot_flux = self.flux.sum()*np.diff(self.ebins)[0]

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

    def generate_photons(self, n_photons, prng=None):
        if prng is None:
            prng = np.random
        randvec = prng.uniform(size=n_photons)
        randvec.sort()
        cumspec = np.cumsum(self.flux)
        cumspec = np.insert(cumspec, 0, 0.0)
        cumspec /= cumspec[-1]
        energies = np.interp(randvec, cumspec, self.ebins)
        return energies

class PowerLawSpectrum(Spectrum):
    def __init__(self, emin, emax, alpha, tot_flux, nbins=10000):
        self.ebins = np.linspace(emin, emax, nbins+1)
        self.emid = 0.5*(self.ebins[1:]+self.ebins[:-1])
        self.tot_flux = tot_flux
        flux = self.emid**alpha
        flux /= flux.sum()
        self.flux = flux*self.tot_flux/np.diff(self.ebins)[0]
        self.nbins = nbins