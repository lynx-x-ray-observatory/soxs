from soxs.spatial import PointSourceModel
from soxs.spectra import Spectrum
from soxs.simput import SimputPhotonList, SimputCatalog
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

    pos1 = PointSourceModel(ra0+0.05, dec0+0.05)
    pos2 = PointSourceModel(ra0-0.05, dec0-0.05)

    sc = SimputCatalog.from_models("pt_src1", (spec, pos1), exp_time, area)
    pl2 = SimputPhotonList.from_models("pt_src2", spec, pos2, exp_time, area)
    sc.append(pl2)
    fns = {0: "pt_src1_phlist.fits", 1: "pt_src2_phlist.fits"}
    sc.write_catalog("pt_src_simput.fits", src_filenames=fns, overwrite=True)

    assert os.path.exists("pt_src1_phlist.fits")
    assert os.path.exists("pt_src2_phlist.fits")
    assert os.path.exists("pt_src_simput.fits")

    f = pyfits.open("pt_src_simput.fits")
    cat = f["SRC_CAT"].data["SPECTRUM"]
    assert cat[0] == "pt_src1_phlist.fits[PHLIST,1]"
    assert cat[1] == "pt_src2_phlist.fits[PHLIST,1]"
    f.close()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

