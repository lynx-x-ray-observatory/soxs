import numpy as np
from soxs.utils import construct_wcs
import astropy.units as u

def generate_radial_events(num_events, func, prng=np.random):
    rbins = np.linspace(0.0, 3000.0, 100000)
    rmid = 0.5*(rbins[1:]+rbins[:-1])
    pdf = func(rmid)*rmid
    pdf /= pdf.sum()
    radius = prng.choice(rmid, size=num_events, p=pdf)
    theta = 2.*np.pi*prng.uniform(size=num_events)
    x = radius*np.cos(theta)
    y = radius*np.sin(theta)
    return x, y

class SpatialModel(object):
    def __init__(self, ra, dec):
        self.ra = u.Quantity(ra, "deg")
        self.dec = u.Quantity(dec, "deg")
        self.num_events = self.ra.size

    def __add__(self, other):
        ra = np.concatenate([self.ra, other.ra])
        dec = np.concatenate([self.dec, other.dec])
        return SpatialModel(ra, dec)

class PointSourceModel(SpatialModel):
    """
    Create positions for a photons emanating from a point source.

    Parameters
    ----------
    ra0 : float
        The RA of the source in degrees.
    dec0 : float
        The Dec of the source in degrees.
    num_events : integer
        The number of events to generate. 
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """

    def __init__(self, ra0, dec0, num_events):
        ra = ra0*np.ones(num_events)
        dec = dec0*np.ones(num_events)
        super(PointSourceModel, self).__init__(ra, dec)

class RadialFunctionModel(SpatialModel):
    """
    Create positions for photons using a generic surface brightness
    profile as a function of radius.

    Parameters
    ----------
    ra0 : float
        The center RA of the source in degrees.
    dec0 : float
        The center Dec of the source in degrees.
    func : function or function-like, something callable.
        A function that takes an array of radii and generates a radial
        surface brightness profile. 
    num_events : integer
        The number of events to generate. 
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, func, num_events, prng=np.random):
        x, y = generate_radial_events(num_events, func,
                                      prng=prng)
        w = construct_wcs(ra0, dec0)
        ra, dec = w.wcs_pix2world(x, y, 1)
        super(RadialFunctionModel, self).__init__(ra, dec)

class RadialArrayModel(RadialFunctionModel):
    """
    Create positions for photons using a table of radii and 
    surface brightness contained in two arrays. 

    Parameters
    ----------
    ra0 : float
        The center RA of the source in degrees.
    dec0 : float
        The center Dec of the source in degrees.
    r : NumPy array
        The array of radii for the profile in arcseconds. 
    S_r: float
        The array of the surface brightness of the profile. 
    num_events : integer
        The number of events to generate. 
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, r, S_r, num_events, prng=np.random):
        func = lambda rr: np.interp(rr, r, S_r, left=0.0, right=0.0)
        super(RadialArrayModel, self).__init__(ra0, dec0, func,
                                               num_events, prng=prng)

class RadialFileModel(RadialArrayModel):
    """
    Create positions for photons using a table of radii and 
    surface brightness contained in a file. 

    Parameters
    ----------
    ra0 : float
        The center RA of the source in degrees.
    dec0 : float
        The center Dec of the source in degrees.
    radfile : string
        The file containing the table of radii and surface 
        brightness. It must be an ASCII table with only two
        columns. 
    num_events : integer
        The number of events to generate. 
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, radfile, num_events, prng=np.random):
        r, S_r = np.loadtxt(radfile, unpack=True)
        super(RadialFileModel, self).__init__(ra0, dec0, r, S_r, 
                                              num_events, prng=prng)

class BetaModel(RadialFunctionModel):
    """
    Create positions for photons with a beta-model shape.

    Parameters
    ----------
    ra0 : float
        The center RA of the beta model in degrees.
    dec0 : float
        The center Dec of the beta model in degrees.
    r_c: float
        The core radius of the profile in arcseconds.
    beta : float
        The "beta" parameter of the profile.
    num_events : integer
        The number of events to generate. 
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, r_c, beta, num_events,
                 prng=np.random):
        func = lambda r: (1.0+(r/r_c)**2)**(-3*beta+0.5)
        super(BetaModel, self).__init__(ra0, dec0, func,
                                        num_events, prng=prng)

class AnnulusModel(RadialFunctionModel):
    """
    Create positions for photons within an annulus shape with uniform surface brightness.

    Parameters
    ----------
    ra0 : float
        The center RA of the annulus in degrees.
    dec0 : float
        The center Dec of the annulus in degrees.
    r_in : float
        The inner radius of the annulus in arcseconds.
    r_out: float
        The outer radius of the annulus in arcseconds.
    num_events : integer
        The number of events to generate. 
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, r_in, r_out, num_events,
                 prng=np.random):
        def func(r):
            f = np.zeros(r.size)
            idxs = np.logical_and(r >= r_in, r < r_out)
            f[idxs] = 1.0
            return f
        super(AnnulusModel, self).__init__(ra0, dec0, func,
                                           num_events, prng=prng)

class FillFOVModel(SpatialModel):
    """
    Create positions for photons which span a field of view.

    Parameters
    ----------
    ra0 : float
        The center RA of the field of view in degrees.
    dec0 : float
        The center Dec of the field of view in degrees.
    fov : float
        The width of the field of view in arcminutes.
    num_events : integer
        The number of events to generate. 
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, fov, num_events, prng=np.random):
        w = construct_wcs(ra0, dec0)
        width = fov*60.0
        x = prng.uniform(low=-0.5*width, high=0.5*width, size=num_events)
        y = prng.uniform(low=-0.5*width, high=0.5*width, size=num_events)
        ra, dec = w.wcs_pix2world(x, y, 1)
        super(FillFOVModel, self).__init__(ra, dec)