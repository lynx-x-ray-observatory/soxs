import h5py
import numpy as np
from scipy.interpolate import interp1d

from soxs.utils import get_data_file


def wabs_cross_section(E):
    emax = np.array(
        [
            0.0,
            0.1,
            0.284,
            0.4,
            0.532,
            0.707,
            0.867,
            1.303,
            1.840,
            2.471,
            3.210,
            4.038,
            7.111,
            8.331,
            10.0,
        ]
    )
    c0 = np.array(
        [
            17.3,
            34.6,
            78.1,
            71.4,
            95.5,
            308.9,
            120.6,
            141.3,
            202.7,
            342.7,
            352.2,
            433.9,
            629.0,
            701.2,
        ]
    )
    c1 = np.array(
        [
            608.1,
            267.9,
            18.8,
            66.8,
            145.8,
            -380.6,
            169.3,
            146.8,
            104.7,
            18.7,
            18.7,
            -2.4,
            30.9,
            25.2,
        ]
    )
    c2 = np.array(
        [
            -2150.0,
            -476.1,
            4.3,
            -51.4,
            -61.1,
            294.0,
            -47.7,
            -31.5,
            -17.0,
            0.0,
            0.0,
            0.75,
            0.0,
            0.0,
        ]
    )
    idxs = np.minimum(np.searchsorted(emax, E) - 1, 13)
    sigma = (c0[idxs] + c1[idxs] * E + c2[idxs] * E * E) * 1.0e-24 / E**3
    return sigma


def get_wabs_absorb(e, nH):
    sigma = wabs_cross_section(e)
    return np.exp(-nH * 1.0e22 * sigma)


_tbabs_emid = None
_tbabs_sigma = None
_tbabs_func = None


def tbabs_cross_section(E, abund_table="angr"):
    global _tbabs_emid
    global _tbabs_sigma
    global _tbabs_func
    if _tbabs_func is None:
        with h5py.File(get_data_file("tbabs_table.h5"), "r") as f:
            _tbabs_sigma = f["cross_section"][abund_table][:]
            nbins = _tbabs_sigma.size
            ebins = np.linspace(f["emin"][()], f["emax"][()], nbins + 1)
        _tbabs_emid = 0.5 * (ebins[1:] + ebins[:-1])
        _tbabs_func = interp1d(
            _tbabs_emid,
            _tbabs_sigma,
            bounds_error=False,
            fill_value=(_tbabs_sigma[0], _tbabs_sigma[-1]),
            assume_sorted=True,
            copy=False,
        )
    return _tbabs_func(E)


def get_tbabs_absorb(e, nH, abund_table="angr"):
    sigma = tbabs_cross_section(e, abund_table=abund_table)
    return np.exp(-nH * 1.0e22 * sigma)
