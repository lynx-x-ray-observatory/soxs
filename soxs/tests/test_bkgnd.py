from soxs.instrument import make_background, AuxiliaryResponseFile, \
    FlatResponse, instrument_simulator, make_background_file
from soxs.background.foreground import hm_astro_bkgnd
from soxs.background.instrument import acisi_particle_bkgnd
from soxs.background.spectra import ConvolvedBackgroundSpectrum
from numpy.random import RandomState
from numpy.testing import assert_allclose
import astropy.io.fits as pyfits
import tempfile
import os
import shutil
import numpy as np

prng = RandomState(24)

def test_uniform_bkgnd_scale():
    hdxi_arf = AuxiliaryResponseFile("xrs_hdxi_3x10.arf")
    flat_arf = FlatResponse(hdxi_arf.elo[0], hdxi_arf.ehi[-1], 1.0, 
                            hdxi_arf.emid.size)
    events, event_params = make_background(50000.0, "hdxi", [30., 45.], 
                                           foreground=True, instr_bkgnd=True,
                                           ptsrc_bkgnd=False, cosmo_bkgnd=False,
                                           prng=prng)
    ncts = np.logical_and(events["energy"] >= 0.7, events["energy"] <= 2.0).sum()
    t_exp = event_params["exposure_time"]
    fov = (event_params["fov"]*60.0)**2
    S = ncts/t_exp/fov
    dS = np.sqrt(ncts)/t_exp/fov
    foreground = ConvolvedBackgroundSpectrum(hm_astro_bkgnd, hdxi_arf)
    instr_bkgnd = ConvolvedBackgroundSpectrum(acisi_particle_bkgnd, flat_arf)
    f_sum = foreground.get_flux_in_band(0.7, 2.0)[0]
    i_sum = instr_bkgnd.get_flux_in_band(0.7, 2.0)[0]
    b_sum = (f_sum+i_sum).to("ph/(arcsec**2*s)").value
    assert np.abs(S-b_sum) < 1.645*dS

def test_add_background():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    prng1 = RandomState(29)
    prng2 = RandomState(29)

    ra0 = 30.0
    dec0 = 45.0
    exp_time = 50000.0

    ra = np.array([])
    dec = np.array([])
    e = np.array([])

    empty_cat = {"ra": [ra], "dec": [dec], "energy": [e],
                 "flux": [0.0], "emin": [0.1], "emax": [10.0],
                 "sources": ["empty"]}

    instrument_simulator(empty_cat, "evt1.fits", exp_time, "hdxi",
                         [ra0, dec0], prng=prng1, clobber=True)

    make_background_file("bkg_evt.fits", exp_time, "hdxi", [ra0, dec0],
                         prng=prng2, clobber=True)

    instrument_simulator(empty_cat, "evt2.fits", exp_time, "hdxi",
                         [ra0, dec0], bkgnd_file="bkg_evt.fits",
                         prng=prng2, clobber=True)

    f1 = pyfits.open("evt1.fits")
    f2 = pyfits.open("evt2.fits")

    for key in ["X", "Y", "ENERGY", "PHA"]:
        assert_allclose(f1["EVENTS"].data[key], f2["EVENTS"].data[key], 
                            )
    f1.close()
    f2.close()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_ptsrc_source_numbers():
    from soxs.background.point_sources import generate_sources
    from soxs.data import cdf_fluxes, cdf_gal, cdf_agn
    prng = RandomState(33)
    fov = 20.0 # arcmin
    exp_time = 500000.0 # seconds
    area = 30000.0 # cm**2
    f_agn = np.zeros((cdf_fluxes.size-1, 100))
    f_gal = np.zeros((cdf_fluxes.size-1, 100))
    for k in range(100):
        agn_sources, gal_sources = generate_sources(exp_time, area, fov, 
                                                    prng=prng)
        agn_fluxes = np.array([agn.flux for agn in agn_sources])
        gal_fluxes = np.array([gal.flux for gal in gal_sources])
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

if __name__ == "__main__":
    test_add_background()
    test_uniform_bkgnd_scale()
    test_ptsrc_source_numbers()