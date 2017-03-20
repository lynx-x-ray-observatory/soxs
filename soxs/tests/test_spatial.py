from soxs.spatial import PointSourceModel, BetaModel, \
    AnnulusModel
from soxs.spectra import ApecGenerator, ConvolvedSpectrum
import numpy as np
import os
import shutil
import tempfile
import astropy.io.fits as pyfits
from soxs.events import write_radial_profile
from soxs.simput import write_photon_list
from soxs.instrument import instrument_simulator, sigma_to_fwhm, \
    AuxiliaryResponseFile
from soxs.instrument_registry import get_instrument_from_registry, \
    add_instrument_to_registry
from numpy.random import RandomState
from sherpa.astro.ui import set_source, freeze, \
    fit, covar, get_covar_results, set_covar_opt, \
    load_data, set_stat, set_method

kT = 6.0
Z = 0.3
redshift = 0.03
norm = 1.0e-3
nH = 0.04
exp_time = 5.0e4
area = 30000.0

prng = RandomState(31)

agen = ApecGenerator(0.05, 12.0, 10000, broadening=True)
spec = agen.get_spectrum(kT, Z, redshift, norm)
spec.apply_foreground_absorption(nH)

ra0 = 30.0
dec0 = 45.0

def test_point_source():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    e = spec.generate_energies(exp_time, area, prng=prng)

    pt_src = PointSourceModel(ra0, dec0, e.size)

    write_photon_list("pt_src", "pt_src", e.flux, pt_src.ra, pt_src.dec,
                      e, overwrite=True)

    inst = get_instrument_from_registry("hdxi")
    inst["name"] = "hdxi_big_psf"
    inst["psf"] = ["gaussian", 5.0]

    add_instrument_to_registry(inst)

    instrument_simulator("pt_src_simput.fits", "pt_src_evt.fits", exp_time,
                         "hdxi_big_psf", [ra0, dec0], ptsrc_bkgnd=False, 
                         instr_bkgnd=False, foreground=False, prng=prng)

    psf_scale = inst["psf"][1]
    dtheta = inst["fov"]*60.0/inst["num_pixels"]

    f = pyfits.open("pt_src_evt.fits")
    x = f["EVENTS"].data["X"]
    y = f["EVENTS"].data["Y"]
    f.close()

    scalex = np.std(x)*sigma_to_fwhm*dtheta
    scaley = np.std(y)*sigma_to_fwhm*dtheta

    assert (scalex - psf_scale)/psf_scale < 0.01
    assert (scaley - psf_scale)/psf_scale < 0.01

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_annulus():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    r_in = 10.0
    r_out = 30.0

    e = spec.generate_energies(exp_time, area, prng=prng)

    ann_src = AnnulusModel(ra0, dec0, r_in, r_out, e.size, prng=prng)

    write_photon_list("ann", "ann", e.flux, ann_src.ra, ann_src.dec,
                      e, overwrite=True)

    instrument_simulator("ann_simput.fits", "ann_evt.fits", exp_time,
                         "hdxi", [ra0, dec0], ptsrc_bkgnd=False, 
                         instr_bkgnd=False, foreground=False, prng=prng)

    inst = get_instrument_from_registry("hdxi")
    arf = AuxiliaryResponseFile(inst["arf"])
    cspec = ConvolvedSpectrum(spec, arf)
    ph_flux = cspec.get_flux_in_band(0.5, 7.0)[0].value
    S0 = ph_flux/(np.pi*(r_out**2-r_in**2))

    write_radial_profile("ann_evt.fits", "ann_evt_profile.fits", [ra0, dec0],
                         1.1*r_in, 0.9*r_out, 100, ctr_type="celestial", 
                         emin=0.5, emax=7.0, overwrite=True)

    load_data(1, "ann_evt_profile.fits", 3, ["RMID","SUR_BRI","SUR_BRI_ERR"])
    set_stat("chi2")
    set_method("levmar")
    set_source("const1d.src")
    src.c0 = 0.8*S0

    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-S0) < res.parmaxes[0]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_beta_model():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    r_c = 20.0
    beta = 1.0

    e = spec.generate_energies(exp_time, area, prng=prng)

    beta_src = BetaModel(ra0, dec0, r_c, beta, e.size, prng=prng)

    write_photon_list("beta", "beta", e.flux, beta_src.ra, beta_src.dec,
                      e, overwrite=True)

    instrument_simulator("beta_simput.fits", "beta_evt.fits", exp_time,
                         "hdxi", [ra0, dec0], ptsrc_bkgnd=False, 
                         instr_bkgnd=False, foreground=False, prng=prng)

    inst = get_instrument_from_registry("hdxi")
    arf = AuxiliaryResponseFile(inst["arf"])
    cspec = ConvolvedSpectrum(spec, arf)
    ph_flux = cspec.get_flux_in_band(0.5, 7.0)[0].value
    S0 = 3.0*ph_flux/(2.0*np.pi*r_c*r_c)

    write_radial_profile("beta_evt.fits", "beta_evt_profile.fits", [ra0, dec0],
                         0.0, 100.0, 200, ctr_type="celestial", emin=0.5, 
                         emax=7.0, overwrite=True)

    load_data(1, "beta_evt_profile.fits", 3, ["RMID","SUR_BRI","SUR_BRI_ERR"])
    set_stat("chi2")
    set_method("levmar")
    set_source("beta1d.src")
    src.beta = 1.0
    src.r0 = 10.0
    src.ampl = 0.8*S0
    freeze(src.xpos)

    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-r_c) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-beta) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-S0) < res.parmaxes[2]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_point_source()
    test_annulus()
    test_beta_model()