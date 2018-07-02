from soxs.spatial import PointSourceModel, BetaModel, \
    AnnulusModel
from soxs.spectra import ApecGenerator
import numpy as np
import os
import shutil
import tempfile
import astropy.io.fits as pyfits
from astropy.units import Quantity
from soxs.events import write_radial_profile, make_exposure_map
from soxs.simput import SimputCatalog
from soxs.instrument import instrument_simulator, sigma_to_fwhm
from soxs.instrument_registry import get_instrument_from_registry, \
    add_instrument_to_registry
from soxs.tests.utils import file_answer_testing

kT = Quantity(6.0, "keV")
Z = 0.3
redshift = 0.03
norm = 1.0e-3
nH = 0.04
exp_time = Quantity(20.0, "ks")
area = Quantity(3.0, "m**2")

prng = 31

agen = ApecGenerator((0.05, "keV"), (12.0, "keV"), 10000, broadening=True)
spec = agen.get_spectrum(kT, Z, redshift, norm)
spec.apply_foreground_absorption(nH)

ra0 = 30.0
dec0 = 45.0


def test_point_source():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    pt_src_pos = PointSourceModel(ra0, dec0)
    sim_cat = SimputCatalog.from_models("pt_src", "pt_src", spec, pt_src_pos,
                                        exp_time, area, prng=prng)
    sim_cat.write_catalog(overwrite=True)

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


def test_annulus(answer_store, answer_dir):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    r_in = 10.0
    r_out = 30.0

    ann_pos = AnnulusModel(ra0, dec0, r_in, r_out)

    sim_cat = SimputCatalog.from_models("ann", "ann", spec, ann_pos,
                                        exp_time, area, prng=prng)
    sim_cat.write_catalog(overwrite=True)

    instrument_simulator("ann_simput.fits", "ann_evt.fits", exp_time,
                         "hdxi", [ra0, dec0], ptsrc_bkgnd=False, 
                         instr_bkgnd=False, foreground=False, prng=prng)

    write_radial_profile("ann_evt.fits", "ann_evt_profile.fits", [ra0, dec0],
                         1.1*r_in, 0.9*r_out, 100, ctr_type="celestial",
                         emin=0.5, emax=7.0, overwrite=True)

    file_answer_testing("EVENTS", "ann_evt.fits", answer_store, answer_dir)
    file_answer_testing("PROFILE", "ann_evt_profile.fits", answer_store, answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_beta_model(answer_store, answer_dir):
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = 32

    r_c = 20.0
    beta = 1.0

    exp_time = Quantity(500.0, "ks")

    beta_src_pos = BetaModel(ra0, dec0, r_c, beta)
    sim_cat = SimputCatalog.from_models("beta", "beta", spec, beta_src_pos,
                                        exp_time, area, prng=prng)
    sim_cat.write_catalog(overwrite=True)

    instrument_simulator("beta_simput.fits", "beta_evt.fits", exp_time,
                         "acisi_cy0", [ra0, dec0], ptsrc_bkgnd=False,
                         instr_bkgnd=False, foreground=False, prng=prng)

    write_radial_profile("beta_evt.fits", "beta_evt_profile.fits", [ra0, dec0],
                         0.0, 100.0, 200, ctr_type="celestial", emin=0.5,
                         emax=7.0, overwrite=True)

    file_answer_testing("EVENTS", "beta_evt.fits", answer_store, answer_dir)
    file_answer_testing("PROFILE", "beta_evt_profile.fits", answer_store, answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_beta_model_flux(answer_store, answer_dir):
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    r_c = 20.0
    beta = 1.0

    prng = 34

    beta_src_pos = BetaModel(ra0, dec0, r_c, beta)
    sim_cat = SimputCatalog.from_models("beta", "beta", spec, beta_src_pos,
                                        exp_time, area, prng=prng)
    sim_cat.write_catalog(overwrite=True)

    instrument_simulator("beta_simput.fits", "beta_flux_evt.fits", exp_time,
                         "acisi_cy0", [ra0, dec0], ptsrc_bkgnd=False,
                         instr_bkgnd=False, foreground=False, 
                         roll_angle=37.0, prng=prng)

    wspec = spec.new_spec_from_band(0.5, 7.0)

    make_exposure_map("beta_flux_evt.fits", "beta_expmap.fits", wspec.emid.value,
                      weights=wspec.flux.value, overwrite=True)

    write_radial_profile("beta_flux_evt.fits", "beta_flux_evt_profile.fits",
                         [ra0, dec0], 0.0, 100.0, 200, ctr_type="celestial",
                         emin=0.5, emax=7.0, expmap_file="beta_expmap.fits",
                         overwrite=True)

    file_answer_testing("EVENTS", "beta_flux_evt.fits", answer_store, answer_dir)
    file_answer_testing("PROFILE", "beta_flux_evt_profile.fits", answer_store, answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
