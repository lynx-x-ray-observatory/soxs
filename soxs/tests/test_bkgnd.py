from soxs.instrument import make_background, \
    instrument_simulator, make_background_file, simulate_spectrum
from soxs.background.diffuse import make_frgnd_spectrum
from soxs.spectra import Spectrum
from soxs.response import AuxiliaryResponseFile, RedistributionMatrixFile
from soxs.utils import soxs_files_path, set_soxs_config
from soxs.tests.utils import spectrum_answer_testing
from numpy.random import RandomState
from numpy.testing import assert_allclose
from astropy.io import fits
import tempfile
import os
import shutil
import numpy as np
import astropy.units as u

acisi_particle_bkgnd = Spectrum.from_file(
    os.path.join(soxs_files_path, "acisi_particle_bkgnd.h5"))


def test_uniform_bkgnd_scale():
    prng = RandomState(25)
    hdxi_arf = AuxiliaryResponseFile("xrs_hdxi_3x10.arf")
    hdxi_rmf = RedistributionMatrixFile("xrs_hdxi.rmf")
    events, event_params = make_background((50, "ks"), "lynx_hdxi", [30., 45.], 
                                           foreground=True, instr_bkgnd=True,
                                           ptsrc_bkgnd=False, prng=prng)
    ch_min = hdxi_rmf.eb_to_ch(0.5)-hdxi_rmf.cmin
    ch_max = hdxi_rmf.eb_to_ch(2.0)-hdxi_rmf.cmin
    ncts = np.logical_and(events[hdxi_rmf.chan_type] >= ch_min, 
                          events[hdxi_rmf.chan_type] <= ch_max).sum()
    t_exp = event_params["exposure_time"]
    fov = (event_params["fov"]*60.0)**2
    S = ncts/t_exp/fov
    dS = np.sqrt(ncts)/t_exp/fov
    foreground = make_frgnd_spectrum(hdxi_arf, hdxi_rmf)
    f_sum = foreground.get_flux_in_band(0.5, 2.0)[0]/u.arcmin**2
    i_sum = acisi_particle_bkgnd.get_flux_in_band(0.5, 2.0)[0]*(u.cm/u.arcmin)**2
    b_sum = (f_sum+i_sum).to_value("ph/(arcsec**2*s)")
    assert np.abs(S-b_sum)/b_sum < 0.02


def test_simulate_bkgnd_spectrum():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng = RandomState(29)

    hdxi_arf = AuxiliaryResponseFile("xrs_hdxi_3x10.arf")
    hdxi_rmf = RedistributionMatrixFile("xrs_hdxi.rmf")

    exp_time = 50000.0
    fov = 3600.0
    simulate_spectrum(None, "lynx_hdxi", exp_time, "test_bkgnd.pha",
                      instr_bkgnd=True, foreground=True, prng=prng,
                      overwrite=True, bkgnd_area=(fov, "arcsec**2"))
    ch_min = hdxi_rmf.eb_to_ch(0.5)-hdxi_rmf.cmin
    ch_max = hdxi_rmf.eb_to_ch(2.0)-hdxi_rmf.cmin
    with fits.open("test_bkgnd.pha") as f:
        ncts = f["SPECTRUM"].data["COUNTS"][ch_min:ch_max].sum()
    S = ncts/exp_time/fov
    dS = np.sqrt(ncts)/exp_time/fov
    foreground = make_frgnd_spectrum(hdxi_arf, hdxi_rmf)
    f_sum = foreground.get_flux_in_band(0.5, 2.0)[0]/u.arcmin**2
    i_sum = acisi_particle_bkgnd.get_flux_in_band(0.5, 2.0)[0]*(u.cm/u.arcmin)**2
    b_sum = (f_sum+i_sum).to_value("ph/(arcsec**2*s)")
    assert np.abs(S-b_sum)/b_sum < 0.02

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_add_background():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng1 = RandomState(29)
    prng2 = RandomState(29)

    ra0 = 30.0
    dec0 = 45.0
    ra1 = 22.0
    dec1 = 22.0
    exp_time = 50000.0

    instrument_simulator(None, "evt1.fits", exp_time, "lynx_hdxi",
                         [ra0, dec0], prng=prng1, overwrite=True)

    make_background_file("bkg_evt.fits", exp_time, "lynx_hdxi", [ra0, dec0],
                         prng=prng2, overwrite=True)

    instrument_simulator(None, "evt2.fits", exp_time, "lynx_hdxi",
                         [ra1, dec1], bkgnd_file="bkg_evt.fits",
                         prng=prng2, overwrite=True)

    with fits.open("evt1.fits") as f1, fits.open("evt2.fits") as f2:
        for key in ["X", "Y", "ENERGY", "PHA"]:
            assert_allclose(f1["EVENTS"].data[key],
                            f2["EVENTS"].data[key], rtol=1.0e-6)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_ptsrc():
    from soxs.background.point_sources import generate_fluxes, \
        make_ptsrc_background
    from soxs.data import cdf_fluxes, cdf_gal, cdf_agn
    from soxs.constants import erg_per_keV
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)
    prng = RandomState(33)
    fov = 20.0
    f_agn = np.zeros((cdf_fluxes.size-1, 100))
    f_gal = np.zeros((cdf_fluxes.size-1, 100))
    for k in range(100):
        agn_fluxes, gal_fluxes = generate_fluxes(fov, prng)
        f_agn[:,k] = np.histogram(agn_fluxes, bins=cdf_fluxes)[0]
        f_gal[:,k] = np.histogram(gal_fluxes, bins=cdf_fluxes)[0]
    mu_agn = np.mean(f_agn, axis=1)
    sigma_agn = np.std(f_agn, axis=1)
    mu_gal = np.mean(f_gal, axis=1)
    sigma_gal = np.std(f_gal, axis=1)
    f_agn0 = np.diff(cdf_agn)*(fov/60.0)**2
    f_gal0 = np.diff(cdf_gal)*(fov/60.0)**2
    err_agn = np.abs(mu_agn-f_agn0)/sigma_agn
    err_gal = np.abs(mu_gal-f_gal0)/sigma_gal
    err_agn[sigma_agn == 0.0] = 0.0
    err_gal[sigma_gal == 0.0] = 0.0
    assert np.all(err_agn < 1.0)
    assert np.all(err_gal < 1.0)
    exp_time = 500000.0 # seconds
    fov = 20.0 # arcmin
    area = 30000.0 # cm**2
    sky_center = [20., 17.]
    prng1 = RandomState(33)
    prng2 = RandomState(33)
    agn_fluxes, gal_fluxes = generate_fluxes(fov, prng2)
    fluxes = np.concatenate([agn_fluxes, gal_fluxes])
    events = make_ptsrc_background(exp_time, fov, sky_center, area=area, 
                                   prng=prng1, nH=-1.0, output_sources="src.dat")
    idxs = np.logical_and(events["energy"] > 0.5, events["energy"] < 2.0)
    E_mean = events["energy"][idxs].mean()*erg_per_keV
    F12 = 0.676e-12  # erg/s/cm**2/deg**2 in 1-2 keV band
    diffuse_flux = 2.0 * F12 * (fov / 60) ** 2
    n1 = (fluxes.sum()+diffuse_flux)*exp_time*area/E_mean
    n2 = idxs.sum()
    dn = np.sqrt(n2)
    assert np.abs(n1-n2) < 1.645*dn
    events2 = make_ptsrc_background(exp_time, fov, sky_center, area=area,
                                    prng=prng1, nH=-1.0, input_sources="src.dat")
    assert_allclose(events["ra"].sum(), events2["ra"].sum(), rtol=1.0e-3)
    assert_allclose(events["dec"].sum(), events2["dec"].sum(), rtol=1.0e-3)
    assert_allclose(events["energy"].sum(), events2["energy"].sum(), rtol=1.0e-3)
    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_change_bkgnd(answer_store, answer_dir):
    set_soxs_config("frgnd_spec_model", "default")
    lem_arf = AuxiliaryResponseFile("lem_030322a.arf")
    lem_rmf = RedistributionMatrixFile("lem_09ev_030322.rmf")
    spectrum_answer_testing(make_frgnd_spectrum(lem_arf, lem_rmf),
                            f"default_frgnd_spectrum.h5", answer_store,
                            answer_dir)
    set_soxs_config("frgnd_spec_model", "halosat")
    spectrum_answer_testing(make_frgnd_spectrum(lem_arf, lem_rmf),
                            f"lem_frgnd_spectrum.h5", answer_store,
                            answer_dir)


if __name__ == "__main__":
    test_add_background()
    test_uniform_bkgnd_scale()
    test_ptsrc()
    test_simulate_bkgnd_spectrum()