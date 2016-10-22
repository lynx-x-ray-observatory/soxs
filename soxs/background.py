import numpy as np
import os
from soxs.spectra import Spectrum
from soxs.utils import check_file_location, \
    construct_wcs, soxs_files_path

class BackgroundSpectrum(Spectrum):
    def __init__(self, filename, fov_scale, bkg_type):
        self.bkg_type = bkg_type
        emid, flux = np.loadtxt(filename, unpack=True)
        de = np.diff(emid)[0]
        ebins = np.append(emid-0.5*de, emid[-1]+0.5*de)
        super(BackgroundSpectrum, self).__init__(ebins, flux)
        self *= 1.0/fov_scale

    def generate_energies(self, t_exp, fov, area=30000.0, prng=None):
        if self.bkg_type == "instrumental":
            A = 1.0
        else:
            A = area
        A *= fov*fov
        return super(BackgroundSpectrum, self).generate_energies(t_exp, A,
                                                                 prng=prng)

acisi_bkgnd_file = os.path.join(soxs_files_path, "acisi_particle_bkgnd.dat")
acisi_particle_bkgnd = BackgroundSpectrum(acisi_bkgnd_file, 282.025, "instrumental")

background_registry = {"acisi": acisi_particle_bkgnd}

def make_astrophysical_background(ra_pnt, dec_pnt, fov, exp_time,
                                  bkgnd_file=None, area=30000.0,
                                  prng=np.random):
    """
    Make events for an astrophysical background, usually for adding to existing
    events.

    Parameters
    ----------
    ra_pnt : float
        The pointing RA of the events.
    dec_pnt : float
        The pointing Dec of the events.
    fov : float
        The field of view on a side, in arcminutes.
    exp_time : float
        The exposure time to use to make the events.
    bkgnd_file : string, optional
        The name of the file to use to make the events containing a spectrum. If
        not supplied, a default astrophysical background supplied with SOXS will
        be used. Default: None
    area : float, optional
        The collecting area used to create the photons, in cm**2. Default: 30000.0
    prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
        A pseudo-random number generator. Typically will only be specified
        if you have a reason to generate the same set of random numbers, such as for a
        test. Default is the :mod:`numpy.random` module.
    """
    events = {}

    if bkgnd_file is None:
        bkgnd_file = check_file_location("hm_cxb_bkgnd.dat", "files")

    bkgnd_spectrum = fov*fov*Spectrum.from_file(bkgnd_file)
    events["energy"] = bkgnd_spectrum.generate_energies(exp_time, area=area, prng=prng)
    n_evt = events["energy"].size

    w = construct_wcs(ra_pnt, dec_pnt)

    width = fov*60.0
    x = prng.uniform(low=-0.5*width, high=0.5*width, size=n_evt)
    y = prng.uniform(low=-0.5*width, high=0.5*width, size=n_evt)
    ra, dec = w.wcs_pix2world(x, y, 1)

    events["ra"] = ra
    events["dec"] = dec

    return events