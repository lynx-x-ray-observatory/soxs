import numpy as np

class Spectrum(object):
    def __init__(self, ebins, flux):
        self.ebins = ebins
        self.emid = 0.5*(ebins[1:]+ebins[:-1])
        self.flux = flux
        self.nbins = len(self.emid)
        self.tot_flux = self.flux.sum()

    @classmethod
    def from_file(cls, filename):
        emid, flux = np.loadtxt(filename, unpack=True)
        de = np.diff(emid)[0]
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