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
    def from_apec(cls, absorb_model, nH, kT, 
                  abund, redshift, norm, broadening=False, 
                  velocity=0.0, emin=0.01, emax=50.0,
                  nbins=10000):
        if broadening:
            model_str = "%s*bapec" % absorb_model
            params = [nH, kT, abund, redshift, velocity, norm]
        else:
            model_str = "%s*apec" % absorb_model
            params = [nH, kT, abund, redshift, norm]
        return cls.from_xspec(model_str, params, emin=emin, emax=emax,
                              nbins=nbins)

    @classmethod
    def from_powerlaw(cls, absorb_model, nH, photon_index,
                      redshift, norm, emin=0.01, emax=50.0,
                      nbins=10000):
        model_str = "%s*powerlaw" % absorb_model
        params = [nH, photon_index, redshift, norm]
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
        if emin is None:
            emin = self.ebins[0]
        if emax is None:
            emax = self.ebins[-1]
        idxs = np.logical_and(self.emid >= emin, self.emid <= emax)
        f = self.flux[idxs].sum()
        self.flux *= new_flux/f
        self.tot_flux = self.flux.sum()

    def generate_photons(self, t_exp, area, prng=None):
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
