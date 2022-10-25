import os
import numpy as np
import shutil
import tempfile
from soxs.response import RedistributionMatrixFile
from soxs.thermal_spectra import ApecGenerator, \
    SpexGenerator, MekalGenerator, CloudyCIEGenerator, \
    IGMGenerator
from soxs.spatial import PointSourceModel
from soxs.simput import SimputCatalog, SimputPhotonList
from soxs.instrument_registry import \
    get_instrument_from_registry
from soxs.instrument import instrument_simulator, \
    simulate_spectrum
from soxs.events import write_spectrum
from numpy.random import RandomState
from numpy.testing import assert_allclose, assert_almost_equal
from soxs.tests.utils import spectrum_answer_testing, \
    file_answer_testing

inst_name = "lynx_lxm"

rmf = RedistributionMatrixFile.from_instrument(inst_name)
agen0 = ApecGenerator(0.01, 10.0, 20000, broadening=True)
agen_var0 = ApecGenerator(0.01, 10.0, 20000, var_elem=["O", "Fe"],
                          broadening=True)
agen_nolines0 = ApecGenerator(0.01, 10.0, 20000, broadening=True,
                              nolines=True)
agen_aspl0 = ApecGenerator(0.01, 10.0, 20000, broadening=True,
                           abund_table="aspl")
agen = ApecGenerator(rmf.elo[0], rmf.ehi[-1], rmf.n_e, broadening=True)
agen_var = ApecGenerator(rmf.elo[0], rmf.ehi[-1], rmf.n_e,
                         var_elem=["O", "Fe"], broadening=True)
agen_nolines = ApecGenerator(rmf.elo[0], rmf.ehi[-1], rmf.n_e,
                             broadening=True, nolines=True)
agen_aspl = ApecGenerator(rmf.elo[0], rmf.ehi[-1], rmf.n_e,
                          broadening=True, abund_table="aspl")
agen_nei = ApecGenerator(rmf.elo[0], rmf.ehi[-1], rmf.n_e,
                         broadening=True, nei=True,
                         var_elem=["O^6","O^3","N^4","Ca^5"])

nH_sim = 0.02
kT_sim = 5.0
abund_sim = 0.4
norm_sim = 1.0e-3
redshift = 0.05
O_sim = 0.4
Fe_sim = 0.4
Ne_sim = 0.4
Si_sim = 0.4
S_sim = 0.4
Mg_sim = 0.4
C_sim = 0.4
N_sim = 0.4
Ca_sim = 0.4

nei_sim = {"O^6": 0.4, "O^3": 0.5, "N^4": 0.7, "Ca^5": 0.9}

exp_time = 5.0e4
area = 40000.0

spec = agen0.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
spec.apply_foreground_absorption(nH_sim)

spec_var = agen_var0.get_spectrum(kT_sim, abund_sim, redshift, norm_sim,
                                  elem_abund={"O": O_sim, "Fe": Fe_sim})
spec_var.apply_foreground_absorption(nH_sim)

spec_nolines = agen_nolines0.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
spec_nolines.apply_foreground_absorption(nH_sim)

spec_aspl = agen_aspl0.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
spec_aspl.apply_foreground_absorption(nH_sim)

spec_nei = agen_nei.get_nei_spectrum(kT_sim, nei_sim, redshift, norm_sim)
spec_nei.apply_foreground_absorption(nH_sim)


def test_thermal(answer_store, answer_dir):

    prng = RandomState(71)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(spec, "thermal_spec.h5", answer_store, answer_dir)

    pt_src_pos = PointSourceModel(30.0, 45.0)
    pt_src = SimputPhotonList.from_models("thermal_model", spec, pt_src_pos,
                                          exp_time, area, prng=prng)
    sim_cat = SimputCatalog.from_source("thermal_model_simput.fits", pt_src,
                                        overwrite=True)

    instrument_simulator("thermal_model_simput.fits", "thermal_model_evt.fits",
                         exp_time, inst_name, [30.0, 45.0], ptsrc_bkgnd=False,
                         foreground=False, instr_bkgnd=False, prng=prng)

    write_spectrum("thermal_model_evt.fits", "thermal_model_evt.pha",
                   overwrite=True)

    file_answer_testing("EVENTS", "thermal_model_evt.fits", answer_store,
                        answer_dir)
    file_answer_testing("SPECTRUM", "thermal_model_evt.pha", answer_store,
                        answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_thermal_from_spectrum(answer_store, answer_dir):

    prng = RandomState(89)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    inst = get_instrument_from_registry(inst_name)

    simulate_spectrum(spec, inst["name"], exp_time,
                      "thermal_model_spec_evt.pha", prng=prng)

    file_answer_testing("SPECTRUM", "thermal_model_spec_evt.pha",
                        answer_store, answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_var_thermal():
    assert_allclose(spec.ebins, spec_var.ebins)
    assert_allclose(spec.flux, spec_var.flux)


def test_nolines_thermal_from_spectrum(answer_store, answer_dir):

    prng = RandomState(101)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    inst = get_instrument_from_registry(inst_name)

    simulate_spectrum(spec_nolines, inst["name"], exp_time,
                      "nolines_thermal_model_evt.pha", prng=prng)

    file_answer_testing("SPECTRUM", "nolines_thermal_model_evt.pha",
                        answer_store, answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_thermal_abund_table(answer_store, answer_dir):

    prng = RandomState(72)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(spec_aspl, "thermal_aspl_spec.h5", answer_store,
                            answer_dir)

    pt_src_pos = PointSourceModel(30.0, 45.0)
    pt_src = SimputPhotonList.from_models("thermal_model_aspl", spec_aspl,
                                          pt_src_pos, exp_time, area, 
                                          prng=prng)
    sim_cat = SimputCatalog.from_source("thermal_model_aspl_simput.fits",
                                        pt_src, overwrite=True)

    instrument_simulator("thermal_model_aspl_simput.fits",
                         "thermal_model_aspl_evt.fits", exp_time, inst_name,
                         [30.0, 45.0], ptsrc_bkgnd=False, foreground=False,
                         instr_bkgnd=False, prng=prng)

    write_spectrum("thermal_model_aspl_evt.fits",
                   "thermal_model_aspl_evt.pha",
                   overwrite=True)

    file_answer_testing("EVENTS", "thermal_model_aspl_evt.fits", answer_store,
                        answer_dir)
    file_answer_testing("SPECTRUM", "thermal_model_aspl_evt.pha", answer_store,
                        answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_thermal_nei(answer_store, answer_dir):

    prng = RandomState(71)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(spec_nei, "thermal_spec_nei.h5", answer_store,
                            answer_dir)

    pt_src_pos = PointSourceModel(30.0, 45.0)
    pt_src = SimputPhotonList.from_models("thermal_model_nei", spec_nei,
                                          pt_src_pos, exp_time, area, prng=prng)
    sim_cat = SimputCatalog.from_source("thermal_model_nei_simput.fits",
                                        pt_src, overwrite=True)

    instrument_simulator("thermal_model_nei_simput.fits",
                         "thermal_model_nei_evt.fits",
                         exp_time, inst_name, [30.0, 45.0], ptsrc_bkgnd=False,
                         foreground=False, instr_bkgnd=False, prng=prng)

    write_spectrum("thermal_model_nei_evt.fits", "thermal_model_nei_evt.pha",
                   overwrite=True)

    file_answer_testing("EVENTS", "thermal_model_nei_evt.fits", answer_store,
                        answer_dir)
    file_answer_testing("SPECTRUM", "thermal_model_nei_evt.pha", answer_store,
                        answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_spex(answer_store, answer_dir):
    spex0 = SpexGenerator(0.01, 10.0, 20000, broadening=True)
    spex_var0 = SpexGenerator(0.01, 10.0, 20000, var_elem=["O", "Fe"],
                              broadening=True)
    specx = spex0.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
    specx.apply_foreground_absorption(nH_sim)

    specx_var = spex_var0.get_spectrum(kT_sim, abund_sim, redshift, norm_sim,
                                       elem_abund={"O": O_sim, "Fe": Fe_sim})
    specx_var.apply_foreground_absorption(nH_sim)

    assert_allclose(specx.ebins, specx_var.ebins)
    assert_allclose(specx.flux, specx_var.flux)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(specx, "spex_spectrum.h5", answer_store,
                            answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_mekal(answer_store, answer_dir):
    mgen = MekalGenerator(0.01, 10.0, 20000)
    mgen_var = MekalGenerator(0.01, 10.0, 20000, var_elem=["O", "Fe"])
    specm = mgen.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)

    specm_var = mgen_var.get_spectrum(kT_sim, abund_sim, redshift, norm_sim,
                                      elem_abund={"O": O_sim, "Fe": Fe_sim})

    assert_allclose(specm.ebins, specm_var.ebins)
    assert_allclose(specm.flux, specm_var.flux)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(specm, "mekal_spectrum.h5", answer_store,
                            answer_dir, rtol=1.0e-5)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_linlog():
    agen0_lin = ApecGenerator(0.01, 10.0, 20000, binscale="linear")
    spec_lin = agen0_lin.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
    agen0_log = ApecGenerator(0.01, 10.0, 20000, binscale="log")
    spec_log = agen0_log.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
    assert_almost_equal(spec_lin.total_flux.value, spec_log.total_flux.value)
    assert_almost_equal(spec_lin.total_energy_flux.value, spec_log.total_energy_flux.value)


def test_cloudy_cie(answer_store, answer_dir):
    cgen = CloudyCIEGenerator(0.5, 10.0, 5000, binscale="log")
    cgen_var1 = CloudyCIEGenerator(0.5, 10.0, 5000, binscale="log",
                                   var_elem=["O", "Ne", "Fe"])
    cgen_var2 = CloudyCIEGenerator(0.5, 10.0, 5000, binscale="log",
                                   var_elem=["O", "Ne", "Fe", "S", "Si", "Mg"])
    cgen_var3 = CloudyCIEGenerator(0.5, 10.0, 5000, binscale="log",
                                   var_elem=["C", "N", "O", "Ne", "Fe", "S", "Si", "Ca", "Mg"])

    cspec = cgen.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
    cspec_var1 = cgen_var1.get_spectrum(kT_sim, abund_sim, redshift, norm_sim,
                                        elem_abund={"O": O_sim, 
                                                    "Ne": Ne_sim,
                                                    "Fe": Fe_sim})
    cspec_var2 = cgen_var2.get_spectrum(kT_sim, abund_sim, redshift, norm_sim,
                                        elem_abund={"O": O_sim,
                                                    "Ne": Ne_sim,
                                                    "Fe": Fe_sim,
                                                    "S": S_sim,
                                                    "Si": Si_sim, 
                                                    "Mg": Mg_sim})
    cspec_var3 = cgen_var3.get_spectrum(kT_sim, abund_sim, redshift, norm_sim,
                                        elem_abund={"C": C_sim,
                                                    "N": N_sim,
                                                    "O": O_sim,
                                                    "Ne": Ne_sim,
                                                    "Fe": Fe_sim,
                                                    "S": S_sim,
                                                    "Si": Si_sim,
                                                    "Ca": Ca_sim,
                                                    "Mg": Mg_sim})

    assert_allclose(cspec.ebins, cspec_var1.ebins)
    assert_allclose(cspec.flux, cspec_var1.flux)

    assert_allclose(cspec.ebins, cspec_var2.ebins)
    assert_allclose(cspec.flux, cspec_var2.flux)

    assert_allclose(cspec.ebins, cspec_var3.ebins)
    assert_allclose(cspec.flux, cspec_var3.flux)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(cspec, "cloudy_spectrum.h5", answer_store,
                            answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_igm(answer_store, answer_dir):
    nH_igm = 1.0e-3
    kT_igm = 0.7

    igen = IGMGenerator(0.2, 5.0, 1000, binscale="log")
    igen_var1 = IGMGenerator(0.2, 5.0, 1000, binscale="log",
                             var_elem=["O", "Ne", "Fe"])
    igen_var2 = IGMGenerator(0.2, 5.0, 1000, binscale="log",
                             var_elem=["O", "Ne", "Fe", "S", "Si", "Mg"])
    igen_var3 = IGMGenerator(0.2, 5.0, 1000, binscale="log",
                             var_elem=["C", "N", "O", "Ne", "Fe", "S", "Si", "Ca", "Mg"])

    ispec = igen.get_spectrum(kT_igm, nH_igm, abund_sim, redshift, norm_sim)
    ispec_var1 = igen_var1.get_spectrum(kT_igm, nH_igm, abund_sim, redshift, norm_sim,
                                        elem_abund={"O": O_sim,
                                                    "Ne": Ne_sim,
                                                    "Fe": Fe_sim})
    ispec_var2 = igen_var2.get_spectrum(kT_igm, nH_igm, abund_sim, redshift, norm_sim,
                                        elem_abund={"O": O_sim,
                                                    "Ne": Ne_sim,
                                                    "Fe": Fe_sim,
                                                    "S": S_sim,
                                                    "Si": Si_sim,
                                                    "Mg": Mg_sim})
    ispec_var3 = igen_var3.get_spectrum(kT_igm, nH_igm, abund_sim, redshift, norm_sim,
                                        elem_abund={"C": C_sim,
                                                    "N": N_sim,
                                                    "O": O_sim,
                                                    "Ne": Ne_sim,
                                                    "Fe": Fe_sim,
                                                    "S": S_sim,
                                                    "Si": Si_sim,
                                                    "Ca": Ca_sim,
                                                    "Mg": Mg_sim})

    assert_allclose(ispec.ebins, ispec_var1.ebins)
    assert_allclose(ispec.flux, ispec_var1.flux, rtol=1.0e-5)

    assert_allclose(ispec.ebins, ispec_var2.ebins)
    assert_allclose(ispec.flux, ispec_var2.flux, rtol=1.0e-5)

    assert_allclose(ispec.ebins, ispec_var3.ebins)
    assert_allclose(ispec.flux, ispec_var3.flux, rtol=1.0e-5)

    sigen = IGMGenerator(0.2, 5.0, 1000, binscale="log", resonant_scattering=True)
    sigen_var1 = IGMGenerator(0.2, 5.0, 1000, binscale="log",
                              var_elem=["O", "Ne", "Fe"],
                              resonant_scattering=True)
    sigen_var2 = IGMGenerator(0.2, 5.0, 1000, binscale="log",
                              var_elem=["O", "Ne", "Fe", "S", "Si", "Mg"],
                              resonant_scattering=True)
    sigen_var3 = IGMGenerator(0.2, 5.0, 1000, binscale="log",
                              var_elem=["C", "N", "O", "Ne", "Fe", "S", "Si", "Ca", "Mg"],
                              resonant_scattering=True)

    sispec = sigen.get_spectrum(kT_igm, nH_igm, abund_sim, redshift, norm_sim)
    sispec_var1 = sigen_var1.get_spectrum(kT_igm, nH_igm, abund_sim, redshift, norm_sim,
                                          elem_abund={"O": O_sim,
                                                      "Ne": Ne_sim,
                                                      "Fe": Fe_sim})
    sispec_var2 = sigen_var2.get_spectrum(kT_igm, nH_igm, abund_sim, redshift, norm_sim,
                                          elem_abund={"O": O_sim,
                                                      "Ne": Ne_sim,
                                                      "Fe": Fe_sim,
                                                      "S": S_sim,
                                                      "Si": Si_sim,
                                                      "Mg": Mg_sim})
    sispec_var3 = sigen_var3.get_spectrum(kT_igm, nH_igm, abund_sim, redshift, norm_sim,
                                          elem_abund={"C": C_sim,
                                                      "N": N_sim,
                                                      "O": O_sim,
                                                      "Ne": Ne_sim,
                                                      "Fe": Fe_sim,
                                                      "S": S_sim,
                                                      "Si": Si_sim,
                                                      "Ca": Ca_sim,
                                                      "Mg": Mg_sim})

    assert_allclose(sispec.ebins, sispec_var1.ebins)
    assert_allclose(sispec.flux, sispec_var1.flux, rtol=1.0e-5)

    assert_allclose(sispec.ebins, sispec_var2.ebins)
    assert_allclose(sispec.flux, sispec_var2.flux, rtol=1.0e-5)

    assert_allclose(sispec.ebins, sispec_var3.ebins)
    assert_allclose(sispec.flux, sispec_var3.flux, rtol=1.0e-5)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(ispec, "igm_spectrum.h5", answer_store,
                            answer_dir, rtol=1.0e-5)

    spectrum_answer_testing(sispec, "igm_scatt_spectrum.h5", answer_store,
                            answer_dir, rtol=1.0e-5)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
