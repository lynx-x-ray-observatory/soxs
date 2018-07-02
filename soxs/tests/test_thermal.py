import numpy as np
import os
import shutil
import tempfile
from soxs.spectra import ApecGenerator, get_wabs_absorb
from soxs.spatial import PointSourceModel
from soxs.simput import SimputCatalog
from soxs.instrument_registry import \
    get_instrument_from_registry
from soxs.instrument import instrument_simulator, \
    RedistributionMatrixFile, simulate_spectrum
from soxs.events import write_spectrum
from numpy.random import RandomState
from numpy.testing import assert_allclose
from soxs.tests.utils import spectrum_answer_testing, \
    file_answer_testing

inst_name = "mucal"

rmf = RedistributionMatrixFile("xrs_%s.rmf" % inst_name)
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


def mymodel(pars, x, xhi=None):
    dx = x[1]-x[0]
    wabs = get_wabs_absorb(x+0.5*dx, pars[0])
    apec = agen.get_spectrum(pars[1], pars[2], pars[3], pars[4])
    eidxs = np.logical_and(rmf.elo >= x[0]-0.5*dx, rmf.elo <= x[-1]+0.5*dx)
    return dx*wabs*apec.flux.value[eidxs]


def mymodel_var(pars, x, xhi=None):
    dx = x[1]-x[0]
    wabs = get_wabs_absorb(x+0.5*dx, pars[0])
    apec = agen_var.get_spectrum(pars[1], pars[2], pars[3], pars[4],
                                 elem_abund={"O": pars[5], "Fe": pars[6]})
    eidxs = np.logical_and(rmf.elo >= x[0]-0.5*dx, rmf.elo <= x[-1]+0.5*dx)
    return dx*wabs*apec.flux.value[eidxs]


def mymodel_nolines(pars, x, xhi=None):
    dx = x[1]-x[0]
    wabs = get_wabs_absorb(x+0.5*dx, pars[0])
    apec = agen_nolines.get_spectrum(pars[1], pars[2], pars[3], pars[4])
    eidxs = np.logical_and(rmf.elo >= x[0]-0.5*dx, rmf.elo <= x[-1]+0.5*dx)
    return dx*wabs*apec.flux.value[eidxs]


def mymodel_aspl(pars, x, xhi=None):
    dx = x[1]-x[0]
    wabs = get_wabs_absorb(x+0.5*dx, pars[0])
    apec_aspl = agen_aspl.get_spectrum(pars[1], pars[2], pars[3], pars[4])
    eidxs = np.logical_and(rmf.elo >= x[0]-0.5*dx, rmf.elo <= x[-1]+0.5*dx)
    return dx*wabs*apec_aspl.flux.value[eidxs]

nH_sim = 0.02
kT_sim = 5.0
abund_sim = 0.4
norm_sim = 1.0e-3
redshift = 0.05
O_sim = 0.4
Fe_sim = 0.4

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


def test_thermal(answer_store):

    prng = RandomState(71)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(spec, "thermal_spec.h5", answer_store)

    pt_src_pos = PointSourceModel(30.0, 45.0)
    sim_cat = SimputCatalog.from_models("thermal_model", "thermal_model", spec,
                                        pt_src_pos, exp_time, area, prng=prng)
    sim_cat.write_catalog(overwrite=True)

    instrument_simulator("thermal_model_simput.fits", "thermal_model_evt.fits", exp_time, 
                         inst_name, [30.0, 45.0], ptsrc_bkgnd=False, foreground=False,
                         instr_bkgnd=False, prng=prng)

    write_spectrum("thermal_model_evt.fits", "thermal_model_evt.pha", overwrite=True)

    file_answer_testing("EVENTS", "thermal_model_evt.fits", answer_store)
    file_answer_testing("SPECTRUM", "thermal_model_evt.pha", answer_store)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_thermal_from_spectrum(answer_store):

    prng = RandomState(89)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    inst = get_instrument_from_registry(inst_name)

    simulate_spectrum(spec, inst["name"], exp_time,
                      "thermal_model_spec_evt.pha", prng=prng)

    file_answer_testing("SPECTRUM", "thermal_model_spec_evt.pha", answer_store)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_var_thermal():
    assert_allclose(spec.ebins, spec_var.ebins)
    assert_allclose(spec.flux, spec_var.flux)


def test_nolines_thermal_from_spectrum(answer_store):

    prng = RandomState(101)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    inst = get_instrument_from_registry(inst_name)

    simulate_spectrum(spec_nolines, inst["name"], exp_time,
                      "nolines_thermal_model_evt.pha", prng=prng)

    file_answer_testing("SPECTRUM", "nolines_thermal_model_evt.pha",
                        answer_store)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_thermal_abund_table(answer_store):

    prng = RandomState(72)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spectrum_answer_testing(spec_aspl, "thermal_aspl_spec.h5", answer_store)

    pt_src_pos = PointSourceModel(30.0, 45.0)
    sim_cat = SimputCatalog.from_models("thermal_model_aspl", "thermal_model_aspl",
                                        spec_aspl, pt_src_pos,
                                        exp_time, area, prng=prng)
    sim_cat.write_catalog(overwrite=True)

    instrument_simulator("thermal_model_aspl_simput.fits",
                         "thermal_model_aspl_evt.fits", exp_time, inst_name,
                         [30.0, 45.0], ptsrc_bkgnd=False, foreground=False,
                         instr_bkgnd=False, prng=prng)

    write_spectrum("thermal_model_aspl_evt.fits",
                   "thermal_model_aspl_evt.pha",
                   overwrite=True)

    file_answer_testing("EVENTS", "thermal_model_aspl_evt.fits", answer_store)
    file_answer_testing("SPECTRUM", "thermal_model_aspl_evt.pha", answer_store)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


if __name__ == "__main__":
    import sys
    answer_store = bool(sys.argv[1])
    test_thermal(answer_store)
    test_thermal_from_spectrum(answer_store)
    test_var_thermal()
    test_nolines_thermal_from_spectrum(answer_store)
    test_thermal_abund_table(answer_store)