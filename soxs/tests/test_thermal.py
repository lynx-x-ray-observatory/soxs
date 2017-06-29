import numpy as np
import os
import shutil
import tempfile
from soxs.spectra import ApecGenerator, get_wabs_absorb
from soxs.spatial import PointSourceModel
from soxs.simput import write_photon_list
from soxs.instrument_registry import \
    get_instrument_from_registry
from soxs.instrument import instrument_simulator, \
    RedistributionMatrixFile, AuxiliaryResponseFile, \
    simulate_spectrum
from soxs.utils import convert_rmf
from soxs.events import write_spectrum
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    covar, get_covar_results, set_covar_opt
from numpy.random import RandomState

inst_name = "mucal"

rmf = RedistributionMatrixFile("xrs_%s.rmf" % inst_name)
agen = ApecGenerator(rmf.elo[0], rmf.ehi[-1], rmf.n_de, broadening=True)
agen_var = ApecGenerator(rmf.elo[0], rmf.ehi[-1], rmf.n_de, 
                         var_elem=["O", "Ca"], broadening=True)
agen_nolines = ApecGenerator(rmf.elo[0], rmf.ehi[-1], rmf.n_de, 
                             broadening=True, nolines=True)

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
                                 elem_abund={"O": pars[5], "Ca": pars[6]})
    eidxs = np.logical_and(rmf.elo >= x[0]-0.5*dx, rmf.elo <= x[-1]+0.5*dx)
    return dx*wabs*apec.flux.value[eidxs]

def mymodel_nolines(pars, x, xhi=None):
    dx = x[1]-x[0]
    wabs = get_wabs_absorb(x+0.5*dx, pars[0])
    apec = agen_nolines.get_spectrum(pars[1], pars[2], pars[3], pars[4])
    eidxs = np.logical_and(rmf.elo >= x[0]-0.5*dx, rmf.elo <= x[-1]+0.5*dx)
    return dx*wabs*apec.flux.value[eidxs]

nH_sim = 0.02
kT_sim = 6.0
abund_sim = 0.4
norm_sim = 1.0e-3
redshift = 0.05
O_sim = 0.3
Ca_sim = 0.5

exp_time = 5.0e4
area = 40000.0

spec = agen.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
spec.apply_foreground_absorption(nH_sim)

spec_var = agen_var.get_spectrum(kT_sim, abund_sim, redshift, norm_sim,
                                 elem_abund={"O": O_sim, "Ca": Ca_sim})
spec_var.apply_foreground_absorption(nH_sim)

spec_nolines = agen_nolines.get_spectrum(kT_sim, abund_sim, redshift, norm_sim)
spec_nolines.apply_foreground_absorption(nH_sim)

def test_thermal():

    prng = RandomState(71)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    e = spec.generate_energies(exp_time, area, prng=prng)

    pt_src = PointSourceModel(30.0, 45.0, e.size)

    write_photon_list("thermal_model", "thermal_model", e.flux, pt_src.ra, pt_src.dec,
                      e, overwrite=True)

    instrument_simulator("thermal_model_simput.fits", "thermal_model_evt.fits", exp_time, 
                         inst_name, [30.0, 45.0], ptsrc_bkgnd=False, foreground=False,
                         instr_bkgnd=False, prng=prng)

    inst = get_instrument_from_registry(inst_name)
    arf = AuxiliaryResponseFile(inst["arf"])
    rmf = RedistributionMatrixFile(inst["rmf"])
    os.system("cp %s ." % arf.filename)
    convert_rmf(rmf.filename)

    write_spectrum("thermal_model_evt.fits", "thermal_model_evt.pha", overwrite=True)

    load_user_model(mymodel, "wapec")
    add_user_pars("wapec", ["nH", "kT", "abund", "redshift", "norm"],
                  [0.01, 4.0, 0.2, redshift, norm_sim*0.8],
                  parmins=[0.0, 0.1, 0.0, -20.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 20.0, 1.0e9],
                  parfrozen=[False, False, False, True, False])

    load_pha("thermal_model_evt.pha")
    set_stat("cstat")
    set_method("simplex")
    set_model("wapec")
    ignore(":0.5, 8.0:")
    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-nH_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-kT_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-abund_sim) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-norm_sim) < res.parmaxes[3]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_thermal_from_spectrum():

    prng = RandomState(65)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    inst = get_instrument_from_registry(inst_name)

    simulate_spectrum(spec, inst["name"], exp_time,
                      "thermal_model_evt.pha", prng=prng)

    arf = AuxiliaryResponseFile(inst["arf"])
    rmf = RedistributionMatrixFile(inst["rmf"])
    os.system("cp %s ." % arf.filename)
    convert_rmf(rmf.filename)

    load_user_model(mymodel, "wapec")
    add_user_pars("wapec", ["nH", "kT", "abund", "redshift", "norm"],
                  [0.01, 4.0, 0.2, redshift, norm_sim*0.8],
                  parmins=[0.0, 0.1, 0.0, -20.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 20.0, 1.0e9],
                  parfrozen=[False, False, False, True, False])

    load_pha("thermal_model_evt.pha")
    set_stat("cstat")
    set_method("simplex")
    set_model("wapec")
    ignore(":0.5, 8.0:")
    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-nH_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-kT_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-abund_sim) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-norm_sim) < res.parmaxes[3]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_var_thermal_from_spectrum():

    prng = RandomState(65)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    inst = get_instrument_from_registry(inst_name)

    simulate_spectrum(spec_var, inst["name"], exp_time,
                      "var_thermal_model_evt.pha", prng=prng)

    arf = AuxiliaryResponseFile(inst["arf"])
    rmf = RedistributionMatrixFile(inst["rmf"])
    os.system("cp %s ." % arf.filename)
    convert_rmf(rmf.filename)

    load_user_model(mymodel_var, "wapec")
    add_user_pars("wapec", ["nH", "kT", "abund", "redshift", "norm", "O", "Ca"],
                  [nH_sim, 4.0, 0.2, redshift, norm_sim*0.8, 0.5, 0.3],
                  parmins=[0.0, 0.1, 0.0, -20.0, 0.0, 0.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 20.0, 1.0e9, 10.0, 10.0],
                  parfrozen=[True, False, False, True, False, False, False])

    load_pha("var_thermal_model_evt.pha")
    set_stat("cstat")
    set_method("simplex")
    set_model("wapec")
    ignore(":0.5, 8.0:")
    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-kT_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-abund_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-norm_sim) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-O_sim) < res.parmaxes[3]
    assert np.abs(res.parvals[4]-Ca_sim) < res.parmaxes[4]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_nolines_thermal_from_spectrum():

    prng = RandomState(65)

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    inst = get_instrument_from_registry(inst_name)

    simulate_spectrum(spec_nolines, inst["name"], exp_time,
                      "nolines_thermal_model_evt.pha", prng=prng)

    arf = AuxiliaryResponseFile(inst["arf"])
    rmf = RedistributionMatrixFile(inst["rmf"])
    os.system("cp %s ." % arf.filename)
    convert_rmf(rmf.filename)

    load_user_model(mymodel_nolines, "wapec")
    add_user_pars("wapec", ["nH", "kT", "abund", "redshift", "norm"],
                  [0.01, 4.0, 0.2, redshift, norm_sim*0.8],
                  parmins=[0.0, 0.1, 0.0, -20.0, 0.0],
                  parmaxs=[10.0, 20.0, 10.0, 20.0, 1.0e9],
                  parfrozen=[False, False, False, True, False])

    load_pha("nolines_thermal_model_evt.pha")
    set_stat("cstat")
    set_method("simplex")
    set_model("wapec")
    ignore(":0.5, 8.0:")
    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-nH_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-kT_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-abund_sim) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-norm_sim) < res.parmaxes[3]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_thermal()
    test_thermal_from_spectrum()
    test_var_thermal_from_spectrum()
    test_nolines_thermal_from_spectrum()