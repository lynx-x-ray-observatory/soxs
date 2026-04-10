import numpy as np
from astropy.units import Quantity
from tqdm import tqdm

from soxs.constants import elem_names
from soxs.utils import mylog


def make_table(cxgen, redshift, ions, elems, kTs, colls):
    nkT = kTs.size
    nc = colls.size
    npts = nkT * nc
    nelem = len(elems)

    # Shift the energy bins to the source frame
    cxgen.model.set_ebins(cxgen.ebins * (1.0 + redshift))

    # collision parameter parsing
    if not hasattr(colls, "unit"):
        colls = Quantity(colls, cxgen.coll_units)
    colls = colls.to_value(cxgen.coll_units)

    if nkT > 1:
        mspec = np.zeros((2, npts, cxgen.nbins))
    else:
        mspec = None
    if nelem > 1:
        vspec = np.zeros((2, nelem, npts, cxgen.nbins))
    else:
        vspec = None

    if mspec is not None:
        abund = np.ones(cxgen.num_elements)
        for e in elems:
            abund[e] = 0.0
        cxgen.model.set_abund(abund * cxgen._atable, elements=cxgen.elements)
        pbar = tqdm(leave=True, total=nkT, desc="Preparing spectrum tables for non-variable elements ")
        for j, kT in enumerate(kTs):
            cxgen.model.set_temperature(kT)
            cxgen.model.set_donorabund(["H", "He"], [1.0, 0.0])
            for k, coll in enumerate(colls):
                mspec[0, k + j * nkT, :] = cxgen.model.calc_spectrum(coll)
            cxgen.model.set_donorabund(["H", "He"], [0.0, 1.0])
            for k, coll in enumerate(colls):
                mspec[1, k + j * nkT, :] = cxgen.model.calc_spectrum(coll)
            pbar.update()
        pbar.close()
        if vspec is not None:
            pbar = tqdm(leave=True, total=nkT, desc="Preparing spectrum tables for variable elements ")
            for j, kT in enumerate(kTs):
                cxgen.model.set_temperature(kT)
                for i, elem in enumerate(elems):
                    abund = np.zeros(cxgen.num_elements)
                    abund[elem] = 1.0
                    cxgen.model.set_abund(abund * cxgen._atable, elements=cxgen.elements)
                    cxgen.model.set_donorabund(["H", "He"], [1.0, 0.0])
                    for k, coll in enumerate(colls):
                        vspec[0, i, k + j * nkT, :] = cxgen.model.calc_spectrum(coll)
                    cxgen.model.set_donorabund(["H", "He"], [0.0, 1.0])
                    for k, coll in enumerate(colls):
                        vspec[1, i, k + j * nkT, :] = cxgen.model.calc_spectrum(coll)
                pbar.update()
            pbar.close()
    else:
        cxgen.model.set_abund(cxgen._atable, elements=cxgen.elements)
        for i, (Z, ion) in enumerate(ions):
            mylog.info("Creating spectrum table for %s %i", elem_names[Z], ion)
            # Set the ionization fraction
            ionfrac = {}
            for e in cxgen.model.elements:
                ionfrac[e] = np.zeros(e + 1)
            ionfrac[Z][ion] = 1.0
            cxgen.model.set_ionfrac(ionfrac)
            cxgen.model.set_donorabund(["H", "He"], [1.0, 0.0])
            for k, coll in enumerate(colls):
                vspec[0, i, k, :] = cxgen.model.calc_spectrum(coll)
            cxgen.model.set_donorabund(["H", "He"], [0.0, 1.0])
            for k, coll in enumerate(colls):
                vspec[1, i, k, :] = cxgen.model.calc_spectrum(coll)

    return mspec, vspec
