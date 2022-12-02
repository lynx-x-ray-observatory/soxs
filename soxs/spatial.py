import numpy as np
from soxs.constants import one_arcsec
from soxs.utils import parse_prng, parse_value, \
    get_rot_mat
import astropy.units as u
from astropy import wcs
from astropy.io import fits


def construct_wcs(ra0, dec0, dtheta=None, nx=None):
    if dtheta is None:
        dtheta = one_arcsec
    if nx is None:
        crpix = [0.0]*2
    else:
        crpix = [0.5*(nx+1)]*2
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = crpix
    w.wcs.cdelt = [-dtheta, dtheta]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2
    return w


def generate_radial_events(num_events, func, prng, ellipticity=1.0):
    rbins = np.linspace(0.0, 3000.0, 100000)
    rmid = 0.5*(rbins[1:]+rbins[:-1])
    pdf = func(rmid)*rmid
    pdf /= pdf.sum()
    radius = prng.choice(rmid, size=num_events, p=pdf)
    theta = 2.*np.pi*prng.uniform(size=num_events)
    x = radius*np.cos(theta)
    y = radius*np.sin(theta)*ellipticity
    return x, y


def rotate_xy(theta, x, y, inverse=False):
    rot_mat = get_rot_mat(theta)
    if inverse:
        rot_mat = np.linalg.inv(rot_mat)
    coords = np.dot(rot_mat, np.array([x, y]))
    return coords


def gen_img_coords(width, nx, theta):
    x, y = np.mgrid[0:nx,0:nx]-0.5*(nx-1)
    x *= width*60.0 / nx
    y *= width*60.0 / nx
    coords = rotate_xy(theta, x.flatten(), y.flatten(), inverse=True)
    return coords


class SpatialModel:
    def __init__(self, ra0, dec0):
        self.ra0 = parse_value(ra0, "deg")
        self.dec0 = parse_value(dec0, "deg")
        self.w = construct_wcs(self.ra0, self.dec0)

    def _generate_coords(self, num_events, prng):
        pass

    def _generate_image(self, width, nx):
        pass

    def generate_image(self, width, nx):
        """
        Generate an ImageHDU from the :class:`~soxs.spatial.SpatialModel`.

        Parameters
        ----------
        width : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The width of the image in arcminutes.
        nx : integer
            The resolution of the image, e.g. the number of pixels
            on a side.
        """
        width = parse_value(width, "arcmin")
        img = self._generate_image(width, nx)
        img /= img.sum()
        dtheta = width/60.0/nx
        w = construct_wcs(0.0, 0.0, dtheta=dtheta, nx=nx)
        imhdu = fits.ImageHDU(data=img, header=w.to_header())
        imhdu.name = "IMAGE"
        imhdu.header["HDUCLASS"] = "HEASARC/SIMPUT"
        imhdu.header["HDUCLAS1"] = "IMAGE"
        imhdu.header["HDUVERS"] = "1.1.0"
        return imhdu

    def generate_coords(self, num_events, prng=None):
        """
        Generate a sample of photon positions from this
        spatial model.

        Parameters
        ----------
        num_events : integer
            The number of events to generate.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only
            be specified if you have a reason to generate the same
            set of random numbers, such as for a test. Default is None,
            which sets the seed based on the system time.
        """
        prng = parse_prng(prng)
        x, y = self._generate_coords(num_events, prng)
        ra, dec = self.w.wcs_pix2world(x, y, 1)
        return u.Quantity(ra, "deg"), u.Quantity(dec, "deg")


class PointSourceModel(SpatialModel):
    """
    A model for positions of photons emanating from 
    a point source.

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The RA of the source in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The Dec of the source in degrees.
    """
    def __init__(self, ra0, dec0):
        super(PointSourceModel, self).__init__(ra0, dec0)

    def _generate_coords(self, num_events, prng):
        return (np.zeros(num_events),)*2
 
    def _generate_image(self, width, nx):
        img = np.zeros((nx,nx))
        img[nx//2, nx//2] = 1.0
        return img


class RadialFunctionModel(SpatialModel):
    """
    A model for positions of photons using a generic 
    surface brightness profile as a function of 
    radius.

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center RA of the source in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center Dec of the source in degrees.
    func : function or function-like, something callable.
        A function that takes an array of radii 
        and generates a radial surface brightness profile. 
    num_events : integer
        The number of events to generate. 
    theta : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The angle through which to rotate the beta model 
        in degrees. Only makes sense if ellipticity is 
        added. Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed 
        as the ratio between the length scales of the x 
        and y coordinates. The value of this parameter 
        will shrink or expand the profile in the direction 
        of the "y" coordinate, so you may need to rotate 
        to get the shape you want. Default: 1.0
    """
    def __init__(self, ra0, dec0, func, theta=0.0, ellipticity=1.0):
        super(RadialFunctionModel, self).__init__(ra0, dec0)
        self.theta = parse_value(theta, "deg")
        self.func = func
        self.ellipticity = ellipticity

    def _generate_coords(self, num_events, prng):
        x, y = generate_radial_events(num_events, self.func, prng,
                                      ellipticity=self.ellipticity)
        coords = rotate_xy(self.theta, x, y)
        return coords[0,:], coords[1,:]

    def _generate_image(self, width, nx):
        coords = gen_img_coords(width, nx, self.theta)
        coords[1,] /= self.ellipticity
        r = np.sqrt((coords**2).sum(axis=0))
        return self.func(r.reshape(nx,nx))


class RadialArrayModel(RadialFunctionModel):
    """
    Create positions for photons using a table of radii and 
    surface brightness contained in two arrays. 

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center RA of the source in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center Dec of the source in degrees.
    r : NumPy array
        The array of radii for the profile in arcseconds. 
    S_r: NumPy array
        The array of the surface brightness of the profile. 
    theta : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The angle through which to rotate the beta model 
        in degrees. Only makes sense if ellipticity is 
        added. Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed 
        as the ratio between the length scales of the x 
        and y coordinates. The value of this parameter will 
        shrink or expand the profile in the direction of the 
        "y" coordinate, so you may need to rotate to get the 
        shape you want. Default: 1.0
    """
    def __init__(self, ra0, dec0, r, S_r, theta=0.0, ellipticity=1.0):
        func = lambda rr: np.interp(rr, r, S_r, left=0.0, right=0.0)
        super(RadialArrayModel, self).__init__(ra0, dec0, func, theta=theta, 
                                               ellipticity=ellipticity)


class RadialFileModel(RadialArrayModel):
    """
    A model for positions of photons using a table of radii and 
    surface brightness contained in a file. 

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center RA of the source in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center Dec of the source in degrees.
    radfile : string
        The file containing the table of radii and surface 
        brightness. It must be an ASCII table with only two
        columns. 
    theta : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The angle through which to rotate the beta model in 
        degrees. Only makes sense if ellipticity is added. 
        Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed 
        as the ratio between the length scales of the x 
        and y coordinates. The value of this parameter will 
        shrink or expand the profile in the direction of the 
        "y" coordinate, so you may need to rotate to get the 
        shape you want. Default: 1.0
    """
    def __init__(self, ra0, dec0, radfile, theta=0.0, ellipticity=1.0):
        r, S_r = np.loadtxt(radfile, unpack=True)
        super(RadialFileModel, self).__init__(ra0, dec0, r, S_r, theta=theta, 
                                              ellipticity=ellipticity)


class BetaModel(RadialFunctionModel):
    """
    A model for positions of photons with a beta-model shape.

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center RA of the beta model in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center Dec of the beta model in degrees.
    r_c: float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The core radius of the profile in arcseconds.
    beta : float
        The "beta" parameter of the profile.
    theta : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The angle through which to rotate the beta model in 
        degrees. Only makes sense if ellipticity is added. 
        Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed 
        as the ratio between the length scales of the x 
        and y coordinates. The value of this parameter will 
        shrink or expand the profile in the direction of the 
        "y" coordinate, so you may need to rotate to get the 
        shape you want. Default: 1.0
    """
    def __init__(self, ra0, dec0, r_c, beta, theta=0.0, 
                 ellipticity=1.0):
        r_c = parse_value(r_c, "arcsec")
        func = lambda r: (1.0+(r/r_c)**2)**(-3*beta+0.5)
        super(BetaModel, self).__init__(ra0, dec0, func, theta=theta, 
                                        ellipticity=ellipticity)


class DoubleBetaModel(RadialFunctionModel):
    """
    A model for positions of photons with a double beta-model shape.

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center RA of the beta model in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center Dec of the beta model in degrees.
    r_c1: float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The inner core radius of the profile in arcseconds.
    beta1 : float
        The inner "beta" parameter of the profile.
    r_c2: float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The outer core radius of the profile in arcseconds.
    beta2 : float
        The outer "beta" parameter of the profile.
    sb_ratio : float
        The ratio of the outer to the inner SB peak value
    theta : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The angle through which to rotate the beta model in 
        degrees. Only makes sense if ellipticity is added. 
        Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed 
        as the ratio between the length scales of the x 
        and y coordinates. The value of this parameter will 
        shrink or expand the profile in the direction of the 
        "y" coordinate, so you may need to rotate to get the 
        shape you want. Default: 1.0
    """
    def __init__(self, ra0, dec0, r_c1, beta1, r_c2, beta2, sb_ratio,
                 theta=0.0, ellipticity=1.0):
        r_c1 = parse_value(r_c1, "arcsec")
        r_c2 = parse_value(r_c2, "arcsec")
        func = lambda r: (1.0+(r/r_c1)**2)**(-3*beta1+0.5) + \
                         sb_ratio*(1.0+(r/r_c2)**2)**(-3*beta2+0.5)
        super(DoubleBetaModel, self).__init__(ra0, dec0, func, theta=theta,
                                              ellipticity=ellipticity)


class AnnulusModel(RadialFunctionModel):
    """
    A model for positions of photons within an annulus shape 
    with uniform surface brightness.

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center RA of the annulus in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center Dec of the annulus in degrees.
    r_in : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The inner radius of the annulus in arcseconds.
    r_out: float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The outer radius of the annulus in arcseconds.
    theta : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The angle through which to rotate the beta model in 
        degrees. Only makes sense if ellipticity is added. 
        Default: 0.0
    ellipticity : float, optional
        The ellipticity of the radial profile, expressed 
        as the ratio between the length scales of the x 
        and y coordinates. The value of this parameter will 
        shrink or expand the profile in the direction of the 
        "y" coordinate, so you may need to rotate to get the 
        shape you want. Default: 1.0
    """
    def __init__(self, ra0, dec0, r_in, r_out, theta=0.0, 
                 ellipticity=1.0):
        r_in = parse_value(r_in, "arcsec")
        r_out = parse_value(r_out, "arcsec")
        def func(r):
            f = np.zeros(r.shape)
            idxs = np.logical_and(r >= r_in, r < r_out)
            f[idxs] = 1.0
            return f
        super(AnnulusModel, self).__init__(ra0, dec0, func, 
                                           theta=theta,
                                           ellipticity=ellipticity)


class RectangleModel(SpatialModel):
    """
    A model for positions of photons within a rectangle 
    or line shape.

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center RA of the rectangle in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center Dec of the rectangle in degrees.
    width : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The width of the rectangle in arcseconds.
    height : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The height of the rectangle in arcseconds.
    theta : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The angle through which to rotate the rectangle 
        in degrees. Default: 0.0
    """
    def __init__(self, ra0, dec0, width, height, theta=0.0):
        super(RectangleModel, self).__init__(ra0, dec0)
        self.width = parse_value(width, "arcsec")
        self.height = parse_value(height, "arcsec")
        self.theta = parse_value(theta, "deg")

    def _generate_image(self, width, nx):
        img = np.zeros((nx, nx))
        c = gen_img_coords(width, nx, self.theta)
        rect = (-0.5*self.width < c[0]) & (c[0] < 0.5*self.width)
        rect &= (-0.5*self.height < c[1]) & (c[1] < 0.5*self.height)
        img[rect.reshape(nx,nx)] = 1.0
        return img

    def _generate_coords(self, num_events, prng):
        x = prng.uniform(low=-0.5*self.width, high=0.5*self.width, size=num_events)
        y = prng.uniform(low=-0.5*self.height, high=0.5*self.height, size=num_events)
        coords = rotate_xy(self.theta, x, y)
        return coords[0,:], coords[1,:]


class FillFOVModel(RectangleModel):
    """
    A model for positions of photons which span a field of view.

    Parameters
    ----------
    ra0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center RA of the field of view in degrees.
    dec0 : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The center Dec of the field of view in degrees.
    fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The width of the field of view in arcminutes.
    """
    def __init__(self, ra0, dec0, fov):
        fov = parse_value(fov, "arcmin")
        width = fov*60.0
        height = fov*60.0
        super(FillFOVModel, self).__init__(ra0, dec0, width, height)
