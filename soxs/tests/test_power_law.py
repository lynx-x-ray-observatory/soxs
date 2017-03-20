import numpy as np
import os
import shutil
import tempfile
from soxs.spectra import Spectrum, get_wabs_absorb
from soxs.spatial import PointSourceModel
from soxs.simput import write_photon_list
from soxs.instrument_registry import \
    get_instrument_from_registry
from soxs.instrument import instrument_simulator, \
    RedistributionMatrixFile, AuxiliaryResponseFile
from soxs.utils import convert_rmf
from soxs.events import write_spectrum
from sherpa.astro.ui import load_user_model, add_user_pars, \
    load_pha, ignore, fit, set_model, set_stat, set_method, \
    covar, get_covar_results, set_covar_opt
from numpy.random import RandomState

prng = RandomState(23)

def mymodel(pars, x, xhi=None):
    dx = x[1]-x[0]
    wabs = get_wabs_absorb(x, pars[0])
    plaw = pars[1]*dx*(x*(1.0+pars[2]))**(-pars[3])
    return wabs*plaw

def test_power_law():
    plaw_fit(1.1)
    plaw_fit(0.8)
    plaw_fit(1.0)

def plaw_fit(alpha_sim):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    nH_sim = 0.02
    norm_sim = 1.0e-4
    redshift = 0.01

    exp_time = 5.0e4
    area = 40000.0
    inst_name = "hdxi"

    spec = Spectrum.from_powerlaw(alpha_sim, redshift, norm_sim)
    spec.apply_foreground_absorption(nH_sim)
    e = spec.generate_energies(exp_time, area, prng=prng)

    pt_src = PointSourceModel(30.0, 45.0, e.size)

    write_photon_list("plaw_model", "plaw_model", e.flux, pt_src.ra, pt_src.dec,
                      e, overwrite=True)

    instrument_simulator("plaw_model_simput.fits", "plaw_model_evt.fits", exp_time, 
                         inst_name, [30.0, 45.0], instr_bkgnd=False, ptsrc_bkgnd=False,
                         foreground=False, prng=prng)

    inst = get_instrument_from_registry(inst_name)
    arf = AuxiliaryResponseFile(inst["arf"])
    rmf = RedistributionMatrixFile(inst["rmf"])
    os.system("cp %s ." % arf.filename)
    convert_rmf(rmf.filename)
    
    write_spectrum("plaw_model_evt.fits", "plaw_model_evt.pha", overwrite=True)

    load_user_model(mymodel, "wplaw")
    add_user_pars("wplaw", ["nH", "norm", "redshift", "alpha"],
                  [0.01, norm_sim*0.8, redshift, 0.9], 
                  parmins=[0.0, 0.0, 0.0, 0.1],
                  parmaxs=[10.0, 1.0e9, 10.0, 10.0],
                  parfrozen=[False, False, True, False])

    load_pha("plaw_model_evt.pha")
    set_stat("cstat")
    set_method("simplex")
    ignore(":0.5, 8.0:")
    set_model("wplaw")
    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-nH_sim) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-norm_sim) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-alpha_sim) < res.parmaxes[2]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_power_law()
