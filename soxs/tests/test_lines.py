import tempfile
import shutil
import os 
from soxs.constants import sigma_to_fwhm, sqrt2pi, \
    hc, ckms
from soxs.events import write_spectrum
from soxs.instrument import instrument_simulator, \
    AuxiliaryResponseFile, RedistributionMatrixFile, \
    simulate_spectrum
from soxs.instrument_registry import \
    get_instrument_from_registry
from soxs.spatial import PointSourceModel
from soxs.simput import SimputCatalog
from soxs.spectra import Spectrum
from soxs.utils import convert_rmf
from numpy.random import RandomState
from sherpa.astro.ui import load_pha, ignore, fit, \
    load_pha, ignore, fit, set_stat, set_method, \
    covar, get_covar_results, set_covar_opt, set_source, \
    load_user_model, add_user_pars, set_model
import numpy as np
from astropy.modeling.functional_models import \
    Gaussian1D

prng = RandomState(69)

def test_emission_line():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    const_flux = 1.0e-4
    line_pos = 5.0
    line_width = 0.02
    line_amp = 1.0e-5

    exp_time = (500.0, "ks")
    area = 40000.0
    inst_name = "mucal"

    spec = Spectrum.from_constant(const_flux, 1.0, 10.0, 20000)
    spec.add_emission_line(line_pos, line_width, line_amp)

    pt_src_pos = PointSourceModel(30.0, 45.0)
    sim_cat = SimputCatalog.from_models("emission_line", "emission_line", spec, pt_src_pos,
                                        exp_time, area, prng=prng)
    sim_cat.write_catalog(overwrite=True)

    instrument_simulator("emission_line_simput.fits", "emission_line_evt.fits", exp_time,
                         inst_name, [30.0, 45.0], instr_bkgnd=False, ptsrc_bkgnd=False,
                         foreground=False, prng=prng)

    inst = get_instrument_from_registry(inst_name)
    arf = AuxiliaryResponseFile(inst["arf"])
    rmf = RedistributionMatrixFile(inst["rmf"])
    os.system("cp %s ." % arf.filename)
    convert_rmf(rmf.filename)

    write_spectrum("emission_line_evt.fits", "emission_line_evt.pha", overwrite=True)

    load_pha("emission_line_evt.pha")
    set_stat("cstat")
    set_method("simplex")
    ignore(":3.0, 7.0:")
    set_source("const1d.src+normgauss1d.ng")
    src.c0 = 1.0e-5
    ng.fwhm = 0.02
    ng.pos = 6.0
    ng.ampl = 1.0e-5
    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-const_flux) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-line_width) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-line_pos) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-line_amp) < res.parmaxes[3]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def mymodel(pars, x, xhi=None):
    dx = x[1]-x[0]
    fwhm = pars[2]*pars[3]/ckms
    sigma = fwhm / sigma_to_fwhm
    B = pars[1] * pars[2] * pars[2]
    B /= hc * sqrt2pi * sigma
    tau = Gaussian1D(B, pars[2], sigma)
    return pars[0]*np.exp(-tau(x+0.5*dx))

def test_absorption_line():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    const_flux = 1.0e-3
    line_pos = 1.0
    line_width = (500.0, "km/s") 
    equiv_width = (3.0e-3, "Angstrom")

    exp_time = (500.0, "ks")
    inst_name = "lynx_gratings"

    spec = Spectrum.from_constant(const_flux, 0.1, 3.0, 200000)
    spec.add_absorption_line(line_pos, line_width, equiv_width)

    simulate_spectrum(spec, inst_name, exp_time, "absorption_line_evt.pha", 
                      overwrite=True)
    inst = get_instrument_from_registry(inst_name)
    arf = AuxiliaryResponseFile(inst["arf"])
    rmf = RedistributionMatrixFile(inst["rmf"])
    os.system("cp %s ." % arf.filename)
    convert_rmf(rmf.filename)

    load_user_model(mymodel, "absorb_line")
    add_user_pars("absorb_line", ["base", "ewidth", "line_center", "sigma"],
                  [3.0e-4, 2.0e-3, 0.8, 400.0],
                  parmins=[1.0e-6, 1.0e-5, 0.1, 10.0],
                  parmaxs=[10.0, 0.1, 10.0, 2000.0],
                  parfrozen=[False, False, False, False])

    load_pha("absorption_line_evt.pha")
    set_stat("cstat")
    set_method("simplex")
    ignore(":0.2, 2.0:")
    set_model("absorb_line")

    fit()
    set_covar_opt("sigma", 1.645)
    covar()
    res = get_covar_results()

    assert np.abs(res.parvals[0]-const_flux) < res.parmaxes[0]
    assert np.abs(res.parvals[1]-equiv_width[0]) < res.parmaxes[1]
    assert np.abs(res.parvals[2]-line_pos) < res.parmaxes[2]
    assert np.abs(res.parvals[3]-line_width[0]) < res.parmaxes[3]

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

if __name__ == "__main__":
    test_emission_line()
    test_absorption_line()