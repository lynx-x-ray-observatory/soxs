import tempfile
import shutil
import os 
from soxs.events import write_spectrum
from soxs.instrument import instrument_simulator, \
    AuxiliaryResponseFile, RedistributionMatrixFile
from soxs.instrument_registry import \
    get_instrument_from_registry
from soxs.spatial import PointSourceModel
from soxs.simput import SimputCatalog
from soxs.spectra import Spectrum
from soxs.utils import convert_rmf
from numpy.random import RandomState
from sherpa.astro.ui import load_pha, ignore, fit, \
    load_pha, ignore, fit, set_stat, set_method, \
    covar, get_covar_results, set_covar_opt, set_source
import numpy as np

prng = RandomState(69)

def test_single_line():
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
    sim_cat = SimputCatalog.from_models("line_model", "line_model", spec, pt_src_pos,
                                        exp_time, area, prng=prng)
    sim_cat.write_catalog(overwrite=True)

    instrument_simulator("line_model_simput.fits", "line_model_evt.fits", exp_time,
                         inst_name, [30.0, 45.0], instr_bkgnd=False, ptsrc_bkgnd=False,
                         foreground=False, prng=prng)

    inst = get_instrument_from_registry(inst_name)
    arf = AuxiliaryResponseFile(inst["arf"])
    rmf = RedistributionMatrixFile(inst["rmf"])
    os.system("cp %s ." % arf.filename)
    convert_rmf(rmf.filename)

    write_spectrum("line_model_evt.fits", "line_model_evt.pha", overwrite=True)

    load_pha("line_model_evt.pha")
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

if __name__ == "__main__":
    test_single_line()