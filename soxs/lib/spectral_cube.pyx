import numpy as np

cimport cython
cimport numpy as np


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def spectral_cube(
    np.ndarray[np.float64_t, ndim=1] x,
    np.ndarray[np.float64_t, ndim=1] y,
    np.ndarray[np.int64_t, ndim=1] c,
    np.ndarray[np.int64_t, ndim=1] cbins,
    int nx, int ny, int nc, int reblock,
):
    cdef double xmin, ymin, dx, dy
    cdef int i, ix, iy, ic, nnx, nny, nnc
    cdef np.ndarray[np.float64_t, ndim=3] data
    cdef int nevent = x.shape[0]

    nnx = int(nx // reblock)
    nny = int(ny // reblock)
    nnc = cbins.size - 1
    data = np.zeros((nnc, nny, nnx), dtype=np.float64)
    dx = float(nx) / nnx
    dy = float(ny) / nny

    cidxs = np.searchsorted(cbins, c) - 1

    for i in range(nevent):
        ix = int((x[i] - 0.5) / dx)
        iy = int((y[i] - 0.5) / dy)
        data[cidxs[i],iy,ix] += 1

    return data
