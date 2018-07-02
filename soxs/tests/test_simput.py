from soxs.simput import write_photon_list
from soxs.spatial import PointSourceModel
from soxs.spectra import Spectrum
import astropy.io.fits as pyfits
import tempfile 
import os
import shutil
from numpy.random import RandomState

prng = RandomState(25)


def test_append():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    exp_time = (50.0, "ks")
    area = (4.0, "m**2")

    ra0 = 30.0
    dec0 = 45.0

    spec = Spectrum.from_powerlaw(1.1, 0.05, 1.0e-4, 0.1, 10.0, 10000)

    e1 = spec.generate_energies(exp_time, area, prng=prng)

    ra1, dec1 = PointSourceModel(ra0+0.05, dec0+0.05).generate_coords(e1.size, prng=prng)

    e2 = spec.generate_energies(exp_time, area, prng=prng)

    ra2, dec2 = PointSourceModel(ra0-0.05, dec0-0.05).generate_coords(e1.size, prng=prng)

    write_photon_list("pt_src", "pt_src1", e1.flux, ra1, dec1, e1, overwrite=True)

    write_photon_list("pt_src", "pt_src2", e2.flux, ra2, dec2, e2, append=True)

    assert os.path.exists("pt_src_simput.fits")
    assert os.path.exists("pt_src1_phlist.fits")
    assert os.path.exists("pt_src2_phlist.fits")

    f = pyfits.open("pt_src_simput.fits")
    cat = f["SRC_CAT"].data["SPECTRUM"]
    assert cat[0] == "pt_src1_phlist.fits[PHLIST,1]"
    assert cat[1] == "pt_src2_phlist.fits[PHLIST,1]"
    f.close()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

