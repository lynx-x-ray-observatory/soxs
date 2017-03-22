import numpy as np
from soxs import write_photon_list
from soxs.constants import keV_per_erg, erg_per_keV
from soxs.spectra import get_wabs_absorb
from soxs.utils import mylog, parse_prng
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.table import Table

# parameters for making event file
# for now, we're leaving these not user-configurable

agn_ind = 1.2  # AGN photon index
gal_ind = 1.2  # galaxy photon index

indices = {"agn": agn_ind, "gal": gal_ind}

fb_emin = 0.5  # keV, low energy bound for the logN-logS flux band
fb_emax = 2.0  # keV, high energy bound for the logN-logS flux band

spec_emin = 0.1  # keV, minimum energy of mock spectrum
spec_emax = 10.0  # keV, max energy of mock spectrum

src_types = ['agn', 'gal']

class Bgsrc:
    def __init__(self, src_type, flux, ind):
        self.src_type = src_type
        self.flux = flux
        self.ind = ind

def get_flux_scale(ind, fb_emin, fb_emax, spec_emin, spec_emax):
    if ind == 1.0:
        f_g = np.log(spec_emax/spec_emin)
    else:
        f_g = (spec_emax**(1.0-ind)-spec_emin**(1.0-ind))/(1.0-ind)
    if ind == 2.0:
        f_E = np.log(fb_emax/fb_emin)
    else:
        f_E = (fb_emax**(2.0-ind)-fb_emin**(2.0-ind))/(2.0-ind)
    fscale = f_g/f_E
    return fscale

def generate_sources(exp_time, area, fov, prng):
    from soxs.data import cdf_fluxes, cdf_gal, cdf_agn

    logf = np.log10(cdf_fluxes)

    n_gal = np.rint(cdf_gal[-1])
    n_agn = np.rint(cdf_agn[-1])
    F_gal = cdf_gal / cdf_gal[-1]
    F_agn = cdf_agn / cdf_agn[-1]
    f_gal = InterpolatedUnivariateSpline(F_gal, logf)
    f_agn = InterpolatedUnivariateSpline(F_agn, logf)

    eph_mean_erg = 1.0*erg_per_keV

    S_min_obs = eph_mean_erg/(exp_time*area)
    mylog.debug("Flux of %g erg/cm^2/s gives roughly "
                "one photon during exposure." % S_min_obs)
    fov_area = fov**2

    n_gal = int(n_gal*fov_area/3600.0)
    n_agn = int(n_agn*fov_area/3600.0)
    mylog.debug("%d AGN, %d galaxies in the FOV." % (n_agn, n_gal))

    randvec1 = prng.uniform(size=n_agn)
    agn_fluxes = 10**f_agn(randvec1)

    randvec2 = prng.uniform(size=n_gal)
    gal_fluxes = 10**f_gal(randvec2)

    agn_sources = [Bgsrc("agn", S, indices["agn"]) for S in agn_fluxes]
    gal_sources = [Bgsrc("gal", S, indices["gal"]) for S in gal_fluxes]

    return agn_sources, gal_sources

def make_ptsrc_background(exp_time, fov, sky_center, nH=0.05, area=40000.0, 
                          output_sources=None, prng=None):
    r"""
    Make a point-source background.

    Parameters
    ----------
    exp_time : float
        The exposure time of the observation in seconds.
    fov : float
        The field of view in arcminutes.
    sky_center : array-like
        The center RA, Dec of the field of view in degrees.
    nH : float, optional
        The hydrogen column in units of 10**22 atoms/cm**2. 
        Default: 0.05
    area : float, optional
        The effective area in cm**2. It must be large enough 
        so that a sufficiently large sample is drawn for the 
        ARF. Default: 40000.
    output_sources : string, optional
        If set to a filename, output the properties of the sources
        within the field of view to a file. Default: None
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    prng = parse_prng(prng)
    agn_sources, gal_sources = generate_sources(exp_time, area, fov, prng)

    sources = agn_sources + gal_sources

    num_sources = len(sources)
    mylog.debug("Generating spectra from %d sources." % num_sources)
    dec_scal = np.fabs(np.cos(sky_center[1]*np.pi/180))
    ra_min = sky_center[0] - fov/(2.0*60.0*dec_scal)
    dec_min = sky_center[1] - fov/(2.0*60.0)
    all_energies = []
    all_ra = []
    all_dec = []

    # This loop is for optimization
    eratio = spec_emax/spec_emin
    oma = {}
    invoma = {}
    fac1 = {}
    fac2 = {}
    for k, ind in indices.items():
        oma[k] = 1.0-ind
        invoma[k] = 1.0/oma[k] if oma[k] != 0.0 else 1.0
        fac1[k] = spec_emin**oma[k]
        fac2[k] = spec_emax**oma[k]-spec_emin**oma[k]

    fluxscale = {}
    for src_type in src_types:
        fluxscale[src_type] = get_flux_scale(indices[src_type], fb_emin,
                                             fb_emax, spec_emin, spec_emax)

    ra0 = prng.uniform(size=num_sources)*fov/(60.0*dec_scal) + ra_min
    dec0 = prng.uniform(size=num_sources)*fov/60.0 + dec_min

    # If requested, output the source properties to a file
    if output_sources is not None:
        t = Table([ra0, dec0, [s.flux for s in sources]],
                  names=('RA', 'Dec', 'flux_0.5_2.0_keV'))
        t["RA"].unit = "deg"
        t["Dec"].unit = "deg"
        t["flux_0.5_2.0_keV"].unit = "erg/(cm**2*s)"
        t.write(output_sources, format='ascii.ecsv', overwrite=True)

    for i, source in enumerate(sources):
        # Using the energy flux, determine the photon flux by simple scaling
        ref_ph_flux = source.flux*fluxscale[source.src_type]*keV_per_erg
        # Now determine the number of photons we will generate
        nph = prng.poisson(ref_ph_flux*exp_time*area)

        if nph > 0:
            # Generate the energies in the source frame
            u = prng.uniform(size=nph)
            if source.ind == 1.0:
                energies = spec_emin*(eratio**u)
            else:
                energies = fac1[source.src_type] + u*fac2[source.src_type]
                energies **= invoma[source.src_type]
            # Assign positions for this source
            ra = ra0[i]*np.ones(nph)
            dec = dec0[i]*np.ones(nph)

            all_energies.append(energies)
            all_ra.append(ra)
            all_dec.append(dec)

    mylog.debug("Finished generating spectra.")

    all_energies = np.concatenate(all_energies)
    all_ra = np.concatenate(all_ra)
    all_dec = np.concatenate(all_dec)

    all_nph = all_energies.size

    # Remove some of the photons due to Galactic foreground absorption.
    # We will throw a lot of stuff away, but this is more general and still
    # faster. 
    if nH is not None:
        absorb = get_wabs_absorb(all_energies, nH)
        randvec = prng.uniform(size=all_energies.size)
        all_energies = all_energies[randvec < absorb]
        all_ra = all_ra[randvec < absorb]
        all_dec = all_dec[randvec < absorb]
        all_nph = all_energies.size
        mylog.debug("%d photons remain after foreground galactic absorption." % all_nph)

    all_flux = np.sum(all_energies)*erg_per_keV/(exp_time*area)

    output_events = {"ra": all_ra, "dec": all_dec, 
                     "energy": all_energies, "flux": all_flux}

    return output_events

def make_point_sources_file(simput_prefix, phlist_prefix, exp_time, fov, 
                            sky_center, nH=0.05, area=40000.0, 
                            prng=None, append=False, overwrite=False,
                            output_sources=None):
    """
    Make a SIMPUT catalog made up of contributions from
    point sources. 

    Parameters
    ----------
    simput_prefix : string
        The filename prefix for the SIMPUT file.
    phlist_prefix : string
        The filename prefix for the photon list file.
    exp_time : float
        The exposure time of the observation in seconds.
    fov : float
        The field of view in arcminutes.
    sky_center : array-like
        The center RA, Dec of the field of view in degrees.
    nH : float, optional
        The hydrogen column in units of 10**22 atoms/cm**2. 
        Default: 0.05
    area : float, optional
        The effective area in cm**2. It must be large enough 
        so that a sufficiently large sample is drawn for the 
        ARF. Default: 40000.
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    append : boolean, optional
        If True, append a new source an existing SIMPUT 
        catalog. Default: False
    overwrite : boolean, optional
        Set to True to overwrite previous files. Default: False
    output_sources : string, optional
        If set to a filename, output the properties of the sources
        within the field of view to a file. Default: None
    """
    events = make_ptsrc_background(exp_time, fov, sky_center, nH=nH, area=area, 
                                   output_sources=output_sources, prng=prng)
    write_photon_list(simput_prefix, phlist_prefix, events["flux"], events["ra"], 
                      events["dec"], events["energy"], append=append, 
                      overwrite=overwrite)
