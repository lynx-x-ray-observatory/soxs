import numpy as np
from xrs_tools.utils import construct_wcs

def generate_radial_events(num_events, func, prng=np.random):
    rbins = np.linspace(0.0, 100.0, 10000)
    pdf = func(rbins)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    randvec = prng.uniform(size=num_events)
    randvec.sort()
    radius = np.interp(randvec, cdf, rbins)
    theta = 2.*np.pi*prng.uniform(size=num_events)
    x = radius*np.cos(theta)
    y = radius*np.sin(theta)
    return x, y

class SpatialModel(object):
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec
        self.num_events = self.ra.size

    def __add__(self, other):
        ra = np.concatenate([self.ra, other.ra])
        dec = np.concatenate([self.dec, other.dec])
        return SpatialModel(ra, dec)

class PointSourceModel(SpatialModel):
    def __init__(self, pt_ra, pt_dec, num_events):
        ra = pt_ra*np.ones(num_events)
        dec = pt_dec*np.ones(num_events)
        super(PointSourceModel, self).__init__(ra, dec)

class RadialFunctionModel(SpatialModel):
    def __init__(self, ra0, dec0, func, num_events, prng=np.random):
        x, y = generate_radial_events(num_events, func,
                                      prng=prng)
        w = construct_wcs(ra0, dec0)
        ra, dec = w.wcs_world2pix(x, y, 1)
        super(RadialFunctionModel, self).__init__(ra, dec)

class RadialArrayModel(RadialFunctionModel):
    def __init__(self, ra0, dec0, r, f_r, num_events, prng=np.random):
        func = lambda rr: np.interp(rr, r, f_r, left=0.0, right=0.0)
        super(RadialArrayModel, self).__init__(ra0, dec0, func,
                                               num_events, prng=prng)

class RadialFileModel(RadialArrayModel):
    def __init__(self, ra0, dec0, radfile, num_events, prng=np.random):
        r, f_r = np.loadtxt(radfile, unpack=True)
        super(RadialFileModel, self).__init__(ra0, dec0, r, f_r, 
                                              num_events, prng=prng)

class BetaModel(RadialFunctionModel):
    def __init__(self, ra0, dec0, beta, r_c, num_events,
                 prng=np.random):
        func = lambda r: (1.0+(r/r_c)**2)**(-3*beta+0.5)
        super(BetaModel, self).__init__(ra0, dec0, func,
                                        num_events, prng=prng)

class AnnulusModel(RadialFunctionModel):
    def __init__(self, ra0, dec0, r_in, r_out, num_events,
                 prng=np.random):
        def func(r):
            f = np.zeros(r.size)
            idxs = np.logical_and(r >= r_in, r < r_out)
            f[idxs] = 1.0
            return f
        super(AnnulusModel, self).__init__(ra0, dec0, func,
                                           num_events, prng=prng)
