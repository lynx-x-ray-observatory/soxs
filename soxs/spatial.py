import numpy as np
from soxs.constants import one_arcsec
import astropy.units as u
import astropy.wcs as pywcs

def construct_wcs(ra0, dec0):
    w = pywcs.WCS(naxis=2)
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [0.0]*2
    w.wcs.cdelt = [-one_arcsec, one_arcsec]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2
    return w

def generate_radial_events(num_events, func, ellipticity=1.0, 
                           prng=np.random):
    rbins = np.linspace(0.0, 3000.0, 100000)
    rmid = 0.5*(rbins[1:]+rbins[:-1])
    pdf = func(rmid)*rmid
    pdf /= pdf.sum()
    radius = prng.choice(rmid, size=num_events, p=pdf)
    theta = 2.*np.pi*prng.uniform(size=num_events)
    x = radius*np.cos(theta)
    y = radius*np.sin(theta)*ellipticity
    return x, y

def rotate_xy(theta, x, y):
    theta_rad = np.deg2rad(theta)
    rot_mat = np.array([[np.cos(theta_rad), -np.sin(theta_rad)],
                        [np.sin(theta_rad), np.cos(theta_rad)]])
    coords = np.dot(rot_mat, np.array([x, y]))
    return coords

class SpatialModel(object):
    def __init__(self, ra, dec, x, y, w):
        self.ra = u.Quantity(ra, "deg")
        self.dec = u.Quantity(dec, "deg")
        self.x = u.Quantity(x, "arcsec")
        self.y = u.Quantity(y, "arcsec")
        self.w = w
        self.num_events = self.ra.size

    def __add__(self, other):
        ra = np.concatenate([self.ra.value, other.ra.value])
        dec = np.concatenate([self.dec.value, other.dec.value])
        x, y = self.w.all_world2pix(ra, dec, 1)
        return SpatialModel(ra, dec, x, y, self.w)

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
        w = construct_wcs(ra0, dec0)
        super(PointSourceModel, self).__init__(ra, dec, 
                                               np.zeros(num_events),
                                               np.zeros(num_events), w)

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
    theta : float, optional
        The angle through which to rotate the beta model in degrees. Only makes
        sense if ellipticity is added. Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed as the ratio between the length
        scales of the x and y coordinates. The value of this parameter will shrink
        or expand the profile in the direction of the "y" coordinate, so you may need 
        to rotate to get the shape you want. Default: 1.0
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, func, num_events, theta=0.0, 
                 ellipticity=1.0, prng=np.random):
        x, y = generate_radial_events(num_events, func,
                                      ellipticity=ellipticity, 
                                      prng=prng)
        w = construct_wcs(ra0, dec0)
        coords = rotate_xy(theta, x, y)
        ra, dec = w.wcs_pix2world(coords[0,:], coords[1,:], 1)
        super(RadialFunctionModel, self).__init__(ra, dec, coords[0,:], 
                                                  coords[1,:], w)

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
    theta : float, optional
        The angle through which to rotate the beta model in degrees. Only makes
        sense if ellipticity is added. Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed as the ratio between the length
        scales of the x and y coordinates. The value of this parameter will shrink
        or expand the profile in the direction of the "y" coordinate, so you may need 
        to rotate to get the shape you want. Default: 1.0
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, r, S_r, num_events, theta=0.0, 
                 ellipticity=1.0, prng=np.random):
        func = lambda rr: np.interp(rr, r, S_r, left=0.0, right=0.0)
        super(RadialArrayModel, self).__init__(ra0, dec0, func, num_events,
                                               theta=theta, ellipticity=ellipticity,
                                               prng=prng)

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
    theta : float, optional
        The angle through which to rotate the beta model in degrees. Only makes
        sense if ellipticity is added. Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed as the ratio between the length
        scales of the x and y coordinates. The value of this parameter will shrink
        or expand the profile in the direction of the "y" coordinate, so you may need 
        to rotate to get the shape you want. Default: 1.0
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, radfile, num_events, theta=0.0, 
                 ellipticity=1.0, prng=np.random):
        r, S_r = np.loadtxt(radfile, unpack=True)
        super(RadialFileModel, self).__init__(ra0, dec0, r, S_r, 
                                              num_events, theta=theta,
                                              ellipticity=ellipticity, prng=prng)

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
    theta : float, optional
        The angle through which to rotate the beta model in degrees. Only makes
        sense if ellipticity is added. Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed as the ratio between the length
        scales of the x and y coordinates. The value of this parameter will shrink
        or expand the profile in the direction of the "y" coordinate, so you may need 
        to rotate to get the shape you want. Default: 1.0
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, r_c, beta, num_events,
                 theta=0.0, ellipticity=1.0, prng=np.random):
        func = lambda r: (1.0+(r/r_c)**2)**(-3*beta+0.5)
        super(BetaModel, self).__init__(ra0, dec0, func,
                                        num_events, theta=theta, 
                                        ellipticity=ellipticity, prng=prng)

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
    theta : float, optional
        The angle through which to rotate the beta model in degrees. Only makes
        sense if ellipticity is added. Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed as the ratio between the length
        scales of the x and y coordinates. The value of this parameter will shrink
        or expand the profile in the direction of the "y" coordinate, so you may need 
        to rotate to get the shape you want. Default: 1.0
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, r_in, r_out, num_events,
                 theta=0.0, ellipticity=1.0, prng=np.random):
        def func(r):
            f = np.zeros(r.size)
            idxs = np.logical_and(r >= r_in, r < r_out)
            f[idxs] = 1.0
            return f
        super(AnnulusModel, self).__init__(ra0, dec0, func, 
                                           num_events, theta=theta,
                                           ellipticity=ellipticity, 
                                           prng=prng)

class RectangleModel(SpatialModel):
    """
    Create positions for photons within a rectangle or line shape.

    Parameters
    ----------
    ra0 : float
        The center RA of the rectangle in degrees.
    dec0 : float
        The center Dec of the rectangle in degrees.
    width : float
        The width of the rectangle in arcseconds.
    height : float
        The height of the rectangle in arcseconds.
    num_events : integer
        The number of events to generate.
    theta : float, optional
        The angle through which to rotate the rectangle in degrees. Default: 0.0
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    def __init__(self, ra0, dec0, width, height, num_events, theta=0.0, prng=np.random):
        w = construct_wcs(ra0, dec0)
        x = prng.uniform(low=-0.5*width, high=0.5*width, size=num_events)
        y = prng.uniform(low=-0.5*height, high=0.5*height, size=num_events)
        coords = rotate_xy(theta, x, y)
        ra, dec = w.wcs_pix2world(coords[0,:], coords[1,:], 1)
        super(RectangleModel, self).__init__(ra, dec, coords[0,:], coords[1,:], w)

class FillFOVModel(RectangleModel):
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
        width = fov*60.0
        height = fov*60.0
        super(FillFOVModel, self).__init__(ra0, dec0, width, height, num_events, prng=prng)
