from soxs.spectra import Spectrum, ConvolvedSpectrum
from soxs.response import AuxiliaryResponseFile
from numpy.testing import assert_allclose, assert_array_equal
import os
import tempfile
import shutil


def test_arithmetic():
    spec1 = Spectrum.from_powerlaw(1.0, 0.05, 1.0e-4, 0.1, 10.0, 10000)
    spec2 = Spectrum.from_powerlaw(2.0, 0.01, 1.0e-3, 0.1, 10.0, 10000)
    spec3 = spec1+spec2
    spec3b = spec1-spec2
    flux3 = spec1.flux+spec2.flux
    flux3b = spec1.flux-spec2.flux

    assert_allclose(spec3.flux, flux3)
    assert_allclose(spec3b.flux, flux3b)

    spec4 = spec3*3.0
    spec5 = 3.0*spec3

    flux4 = spec3.flux*3.0

    assert_allclose(spec4.flux, spec5.flux)
    assert_allclose(spec4.flux, flux4)

    spec6 = spec3/2.5
    flux6 = spec3.flux/2.5

    assert_allclose(spec6.flux, flux6)

    spec7 = Spectrum.from_constant(1.0e-4, 0.1, 10.0, 10000)
    spec8 = spec1+spec7

    assert_allclose(spec8.flux.value, spec1.flux.value+1.0e-4)

    spec8 += spec1

    assert_allclose(spec8.flux.value, (2.0*spec1+spec7).flux.value)


def test_read_write():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    spec1 = Spectrum.from_powerlaw(1.0, 0.05, 1.0e-4, 0.1, 10.0, 10000)
    spec1.write_file("test_spec.dat", overwrite=True)
    spec1.write_file("test_spec.h5", overwrite=True)
    spec1.write_file("test_spec.fits", overwrite=True)

    spec2 = Spectrum.from_file("test_spec.dat")
    spec3 = Spectrum.from_file("test_spec.h5")
    spec4 = Spectrum.from_file("test_spec.fits")
 
    assert_allclose(spec1.flux, spec2.flux)
    assert_allclose(spec1.emid, spec2.emid)
    assert_allclose(spec1.ebins, spec2.ebins)
    assert_allclose(spec1.cumspec, spec2.cumspec)

    assert_allclose(spec1.flux, spec3.flux)
    assert_allclose(spec1.emid, spec3.emid)
    assert_allclose(spec1.ebins, spec3.ebins)
    assert_allclose(spec1.cumspec, spec3.cumspec)

    assert_allclose(spec1.flux, spec4.flux)
    assert_allclose(spec1.emid, spec4.emid)
    assert_allclose(spec1.ebins, spec4.ebins)
    assert_allclose(spec1.cumspec, spec4.cumspec)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_rescale_flux():
    spec = Spectrum.from_powerlaw(2.0, 0.01, 1.0, 0.1, 10.0, 10000)

    spec.rescale_flux(1.0e-4, emin=0.5, emax=7.0, flux_type="photons")
    f = spec.get_flux_in_band(0.5, 7.0)[0]
    assert_allclose(1.0e-4, f.value)

    spec.rescale_flux(1.0e-12, emin=0.4, emax=1.0, flux_type="energy")
    f = spec.get_flux_in_band(0.4, 1.0)[1]
    assert_allclose(1.0e-12, f.value)


def test_convolved_spectra():
    arf = AuxiliaryResponseFile("xrs_hdxi_3x10.arf")
    spec1 = Spectrum.from_powerlaw(2.0, 0.01, 1.0, 0.1, 10.0, 1000)
    cspec1 = ConvolvedSpectrum.convolve(spec1, arf)
    cspec2 = spec1*arf
    spec2 = cspec1.deconvolve()
    assert_array_equal(cspec1.ebins.value, cspec2.ebins.value)
    assert_array_equal(spec1.ebins.value, spec2.ebins.value)
    assert_array_equal(cspec1.flux.value, cspec2.flux.value)
    assert_allclose(spec1.flux.value, spec2.flux.value)