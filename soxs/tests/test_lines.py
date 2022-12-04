import os
import shutil
import tempfile

from soxs.events import write_spectrum
from soxs.instrument import instrument_simulator, simulate_spectrum
from soxs.simput import SimputCatalog, SimputPhotonList
from soxs.spatial import PointSourceModel
from soxs.spectra import Spectrum
from soxs.tests.utils import file_answer_testing, spectrum_answer_testing


def test_emission_line(answer_store, answer_dir):
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    const_flux = 1.0e-4
    line_pos = 5.0
    line_width = 0.02
    line_amp = 1.0e-5

    exp_time = (100.0, "ks")
    area = 30000.0
    inst_name = "lynx_lxm"

    spec = Spectrum.from_constant(const_flux, 1.0, 10.0, 20000)
    spec.add_emission_line(line_pos, line_width, line_amp)

    spectrum_answer_testing(spec, "emission_line_test.h5", answer_store, answer_dir)

    pt_src_pos = PointSourceModel(30.0, 45.0)
    pt_src = SimputPhotonList.from_models(
        "emission_line", spec, pt_src_pos, exp_time, area, prng=69
    )
    SimputCatalog.from_source("emission_line_simput.fits", pt_src, overwrite=True)

    instrument_simulator(
        "emission_line_simput.fits",
        "emission_line_evt.fits",
        exp_time,
        inst_name,
        [30.0, 45.0],
        instr_bkgnd=False,
        ptsrc_bkgnd=False,
        foreground=False,
        prng=69,
    )

    write_spectrum("emission_line_evt.fits", "emission_line_evt.pha", overwrite=True)

    file_answer_testing("EVENTS", "emission_line_evt.fits", answer_store, answer_dir)
    file_answer_testing("SPECTRUM", "emission_line_evt.pha", answer_store, answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_absorption_line(answer_store, answer_dir):
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    const_flux = 1.0e-3
    line_pos = 1.0
    line_width = 0.02
    line_amp = 1.0e-5

    exp_time = (100.0, "ks")
    inst_name = "lynx_xgs"

    spec = Spectrum.from_constant(const_flux, 0.1, 3.0, 100000)
    spec.add_absorption_line(line_pos, line_width, line_amp)

    spectrum_answer_testing(spec, "absorption_line_test.h5", answer_store, answer_dir)

    simulate_spectrum(
        spec, inst_name, exp_time, "absorption_line_evt.pha", overwrite=True, prng=69
    )

    file_answer_testing("SPECTRUM", "absorption_line_evt.pha", answer_store, answer_dir)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
