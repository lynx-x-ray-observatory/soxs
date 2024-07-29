import os
import shutil
import tempfile

import pytest

from soxs.spectra.charge_exchange import ACX2Generator, OneACX2Generator
from soxs.tests.utils import spectrum_answer_testing

acx2 = pytest.importorskip("acx2")

emin = 0.2
emax = 2.0
nbins = 6000
kT = 1.0
collnpar1 = 1.0
collnpar2 = 100.0
abund = 1.0
He_frac = 0.09
C = 0.4
Ne = 0.5
Ca = 0.3
redshift = 0.05
norm = 1.0


def test_vacx2_ctype1(answer_store):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    cxgen0 = ACX2Generator(emin, emax, nbins, collntype=1)

    spec0 = cxgen0.get_spectrum(kT, collnpar1, abund, He_frac, redshift, norm)

    spectrum_answer_testing(spec0, "cx_spec0.h5", answer_store)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_vacx2_ctype2(answer_store):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    cxgen1 = ACX2Generator(emin, emax, nbins, collntype=2)

    spec1 = cxgen1.get_spectrum(kT, collnpar2, abund, He_frac, redshift, norm)

    spectrum_answer_testing(spec1, "cx_spec1.h5", answer_store)

    cxgen2 = ACX2Generator(emin, emax, nbins, collntype=2, var_elem=["C", "Ne", "Ca"])

    spec2 = cxgen2.get_spectrum(
        kT,
        collnpar2,
        abund,
        He_frac,
        redshift,
        norm,
        elem_abund={"C": C, "Ne": Ne, "Ca": Ca},
    )

    spectrum_answer_testing(spec2, "cx_spec2.h5", answer_store)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


def test_oneacx2(answer_store):

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    cxgen4 = OneACX2Generator(emin, emax, nbins, collntype=2)

    spec4 = cxgen4.get_spectrum("Si", 13, collnpar2, He_frac, redshift, norm)

    spectrum_answer_testing(spec4, "cx_spec4.h5", answer_store)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
