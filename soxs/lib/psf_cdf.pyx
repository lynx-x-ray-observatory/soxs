import numpy as np

cimport cython
cimport numpy as np
from libc.math cimport fabs


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def score_psf(np.ndarray[np.float64_t, ndim=1] e_bins,
              np.ndarray[np.float64_t, ndim=1] r_bins,
              np.ndarray[np.float64_t, ndim=1] e,
              np.ndarray[np.float64_t, ndim=1] r):

    # Inspiration here taken from SIMX

    cdef int i, j
    cdef int n_bins = r_bins.size
    cdef int n_evt = r.size
    cdef np.ndarray[np.int64_t, ndim=1] idx_score
    cdef double bestscore = 1e20
    cdef double score

    idx_score = np.zeros(n_evt, dtype='int64')

    for i in range(n_evt):
        for j in range(n_bins):
            score = 100*fabs(e[i]-e_bins[j])+fabs(r[i]-r_bins[j])
            if score < bestscore:
                bestscore = score
                idx_score[i] = j

    return idx_score

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def eef_cdf(np.ndarray[np.int64_t, ndim=1] idx_score,
            np.ndarray[np.float64_t, ndim=1] randvec,
            np.ndarray[np.float64_t, ndim=2] radii,
            np.ndarray[np.float64_t, ndim=2] eef):

    # Binary search portion inspired by:
    # https://eric-bunch.github.io/blog/cython-examples-sampling-and-lda
    cdef int i, j, k
    cdef int kmid, kmax
    cdef double u, de_inv, ep, em
    cdef int n_evt = idx_score.size
    cdef int n_bins = radii.shape[1]
    cdef np.ndarray[np.float64_t, ndim=1] r

    r = np.zeros(n_evt)

    for i in range(n_evt):
        j = idx_score[i]
        u = randvec[i]
        k = 0
        kmax = n_bins
        kmax = n_bins
        while k < kmax:
            kmid = k + ((kmax - k) >> 1)
            if u > eef[j, kmid]:
                k = kmid + 1
            else:
                kmax = kmid
        de_inv = 1.0 / (eef[j, k] - eef[j, k-1])
        ep = (u - eef[j, k-1]) * de_inv
        em = (eef[j, k] - u) * de_inv
        r[i] = radii[j, k-1] * em + radii[j, k] * ep

    return r
