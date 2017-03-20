from soxs.spectra import Spectrum
from numpy.testing import assert_allclose
import os
import tempfile
import shutil

def test_arithmetic():
    spec1 = Spectrum.from_powerlaw(1.0, 0.05, 1.0e-4)
    spec2 = Spectrum.from_powerlaw(2.0, 0.01, 1.0e-3)
    spec3 = spec1+spec2
    flux3 = spec1.flux+spec2.flux

    assert_allclose(spec3.flux, flux3)

    spec4 = spec3*3.0
    spec5 = 3.0*spec3

    flux4 = spec3.flux*3.0

    assert_allclose(spec4.flux, spec5.flux)
    assert_allclose(spec4.flux, flux4)

    spec6 = spec3/2.5
    flux6 = spec3.flux/2.5

    assert_allclose(spec6.flux, flux6)

    spec7 = Spectrum.from_constant(1.0e-4)
    spec8 = spec1+spec7

    assert_allclose(spec8.flux.value, spec1.flux.value+1.0e-4)

def test_read_write():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spec1 = Spectrum.from_powerlaw(1.0, 0.05, 1.0e-4)
    spec1.write_file("test_spec.dat", overwrite=True)
    spec2 = Spectrum.from_file("test_spec.dat")

    assert_allclose(spec1.flux, spec2.flux)
    assert_allclose(spec1.emid, spec2.emid)
    assert_allclose(spec1.ebins, spec2.ebins)
    assert_allclose(spec1.cumspec, spec2.cumspec)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_rescale_flux():
    spec = Spectrum.from_powerlaw(2.0, 0.01, 1.0)

    spec.rescale_flux(1.0e-4, emin=0.5, emax=7.0, flux_type="photons")
    f = spec.get_flux_in_band(0.5, 7.0)[0]
    assert_allclose(1.0e-4, f.value)

    spec.rescale_flux(1.0e-12, emin=0.4, emax=1.0, flux_type="energy")
    f = spec.get_flux_in_band(0.4, 1.0)[1]
    assert_allclose(1.0e-12, f.value)

