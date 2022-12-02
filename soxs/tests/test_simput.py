import os
import shutil
import tempfile

from astropy.io import fits
from numpy.random import RandomState

from soxs.simput import SimputCatalog, SimputPhotonList
from soxs.spatial import PointSourceModel
from soxs.spectra import Spectrum

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

    pos1 = PointSourceModel(ra0 + 0.05, dec0 + 0.05)
    pos2 = PointSourceModel(ra0 - 0.05, dec0 - 0.05)

    pl1 = SimputPhotonList.from_models("pt_src1", spec, pos1, exp_time, area)
    pl2 = SimputPhotonList.from_models("pt_src2", spec, pos2, exp_time, area)
    sc = SimputCatalog.from_source(
        "pt_src_simput.fits", pl1, src_filename="pt_src1_phlist.fits", overwrite=True
    )

    sc.append(pl2, src_filename="pt_src2_phlist.fits", overwrite=True)

    assert os.path.exists("pt_src1_phlist.fits")
    assert os.path.exists("pt_src2_phlist.fits")
    assert os.path.exists("pt_src_simput.fits")

    with fits.open("pt_src_simput.fits") as f:
        cat = f["SRC_CAT"].data["SPECTRUM"].copy()
        if os.name == "nt":
            mypwd = ".\\"
        else:
            mypwd = "./"
        assert cat[0] == f"{mypwd}pt_src1_phlist.fits[PHLIST,1]"
        assert cat[1] == f"{mypwd}pt_src2_phlist.fits[PHLIST,1]"

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
