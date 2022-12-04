import os
import shutil
import tempfile

import numpy as np
from astropy.io import fits
from astropy.units import Quantity

from soxs.constants import sigma_to_fwhm
from soxs.events import make_exposure_map, write_radial_profile
from soxs.instrument import instrument_simulator
from soxs.instrument_registry import (
    add_instrument_to_registry,
    get_instrument_from_registry,
)
from soxs.simput import SimputCatalog, SimputPhotonList, SimputSpectrum
from soxs.spatial import AnnulusModel, BetaModel, DoubleBetaModel
from soxs.tests.utils import file_answer_testing
from soxs.thermal_spectra import ApecGenerator

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

    pt_src = SimputSpectrum.from_spectrum("pt_src", spec, ra0, dec0)
    SimputCatalog.from_source("pt_src_simput.fits", pt_src, overwrite=True)

    inst = get_instrument_from_registry("lynx_hdxi")
    inst["name"] = "hdxi_big_psf"
    inst["psf"] = ["gaussian", 5.0]

    add_instrument_to_registry(inst)

    instrument_simulator(
        "pt_src_simput.fits",
        "pt_src_evt.fits",
        exp_time,
        "hdxi_big_psf",
        [ra0, dec0],
        ptsrc_bkgnd=False,
        instr_bkgnd=False,
        foreground=False,
        prng=prng,
    )

    psf_scale = inst["psf"][1]
    dtheta = inst["fov"] * 60.0 / inst["num_pixels"]

    with fits.open("pt_src_evt.fits") as f:
        x = f["EVENTS"].data["X"].copy()
        y = f["EVENTS"].data["Y"].copy()

    scalex = np.std(x) * sigma_to_fwhm * dtheta
    scaley = np.std(y) * sigma_to_fwhm * dtheta

    assert (scalex - psf_scale) / psf_scale < 0.03
    assert (scaley - psf_scale) / psf_scale < 0.03

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_annulus(answer_store, answer_dir):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    r_in = 10.0
    r_out = 30.0

    ann_pos = AnnulusModel(ra0, dec0, r_in, r_out)

    ann_src = SimputPhotonList.from_models(
        "ann", spec, ann_pos, exp_time, area, prng=prng
    )
    SimputCatalog.from_source("ann_simput.fits", ann_src, overwrite=True)

    instrument_simulator(
        "ann_simput.fits",
        "ann_evt.fits",
        exp_time,
        "lynx_hdxi",
        [ra0, dec0],
        ptsrc_bkgnd=False,
        instr_bkgnd=False,
        foreground=False,
        prng=prng,
    )

    write_radial_profile(
        "ann_evt.fits",
        "ann_evt_profile.fits",
        [ra0, dec0],
        1.1 * r_in,
        0.9 * r_out,
        100,
        ctr_type="celestial",
        emin=0.5,
        emax=7.0,
        overwrite=True,
    )

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
    beta_src = SimputPhotonList.from_models(
        "beta", spec, beta_src_pos, exp_time, area, prng=prng
    )
    SimputCatalog.from_source("beta_simput.fits", beta_src, overwrite=True)

    instrument_simulator(
        "beta_simput.fits",
        "beta_evt.fits",
        exp_time,
        "chandra_acisi_cy0",
        [ra0, dec0],
        ptsrc_bkgnd=False,
        instr_bkgnd=False,
        foreground=False,
        prng=prng,
    )

    write_radial_profile(
        "beta_evt.fits",
        "beta_evt_profile.fits",
        [ra0, dec0],
        0.0,
        100.0,
        200,
        ctr_type="celestial",
        emin=0.5,
        emax=7.0,
        overwrite=True,
    )

    file_answer_testing("EVENTS", "beta_evt.fits", answer_store, answer_dir)
    file_answer_testing("PROFILE", "beta_evt_profile.fits", answer_store, answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_double_beta_model(answer_store, answer_dir):
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = 32

    r_c1 = 20.0
    beta1 = 1.0
    r_c2 = 100.0
    beta2 = 2.0 / 3.0
    sb_ratio = 0.5

    exp_time = Quantity(500.0, "ks")

    double_beta_src_pos = DoubleBetaModel(ra0, dec0, r_c1, beta1, r_c2, beta2, sb_ratio)
    double_beta_src = SimputPhotonList.from_models(
        "double_beta", spec, double_beta_src_pos, exp_time, area, prng=prng
    )
    SimputCatalog.from_source(
        "double_beta_simput.fits", double_beta_src, overwrite=True
    )

    instrument_simulator(
        "double_beta_simput.fits",
        "double_beta_evt.fits",
        exp_time,
        "chandra_acisi_cy0",
        [ra0, dec0],
        ptsrc_bkgnd=False,
        instr_bkgnd=False,
        foreground=False,
        prng=prng,
    )

    write_radial_profile(
        "double_beta_evt.fits",
        "double_beta_evt_profile.fits",
        [ra0, dec0],
        0.0,
        200.0,
        200,
        ctr_type="celestial",
        emin=0.5,
        emax=7.0,
        overwrite=True,
    )

    file_answer_testing("EVENTS", "double_beta_evt.fits", answer_store, answer_dir)
    file_answer_testing(
        "PROFILE", "double_beta_evt_profile.fits", answer_store, answer_dir
    )

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
    beta_src = SimputPhotonList.from_models(
        "beta", spec, beta_src_pos, exp_time, area, prng=prng
    )
    SimputCatalog.from_source("beta_simput.fits", beta_src, overwrite=True)

    instrument_simulator(
        "beta_simput.fits",
        "beta_flux_evt.fits",
        exp_time,
        "chandra_acisi_cy0",
        [ra0, dec0],
        ptsrc_bkgnd=False,
        instr_bkgnd=False,
        foreground=False,
        roll_angle=37.0,
        prng=prng,
    )

    wspec = spec.new_spec_from_band(0.5, 7.0)

    make_exposure_map(
        "beta_flux_evt.fits",
        "beta_expmap.fits",
        wspec.emid.value,
        weights=wspec.flux.value,
        overwrite=True,
    )

    write_radial_profile(
        "beta_flux_evt.fits",
        "beta_flux_evt_profile.fits",
        [ra0, dec0],
        0.0,
        100.0,
        200,
        ctr_type="celestial",
        emin=0.5,
        emax=7.0,
        expmap_file="beta_expmap.fits",
        overwrite=True,
    )

    file_answer_testing("EVENTS", "beta_flux_evt.fits", answer_store, answer_dir)
    file_answer_testing(
        "PROFILE", "beta_flux_evt_profile.fits", answer_store, answer_dir
    )

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
