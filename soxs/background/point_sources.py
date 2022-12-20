import numpy as np
from astropy.io import ascii
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import erf

from soxs.constants import erg_per_keV, keV_per_erg
from soxs.simput import SimputCatalog, SimputPhotonList
from soxs.spatial import construct_wcs
from soxs.spectra import get_tbabs_absorb, get_wabs_absorb
from soxs.utils import mylog, parse_prng, parse_value, soxs_cfg

# Function for computing spectral index of AGN sources
# a fit to the data from Figure 13a of Hickox & Markevitch 2006
# http://adsabs.harvard.edu/abs/2006ApJ...645...95H

# Parameters

aa = -14.0
bb = 0.5
cc = 0.5
dd = 1.8


# Here x = log10(flux)
def get_agn_index(x):
    y = (x - aa) / bb
    return cc * erf(y) + dd


# Index for galaxies
gal_index = 2.0

fb_emin = 0.5  # keV, low energy bound for the logN-logS flux band
fb_emax = 2.0  # keV, high energy bound for the logN-logS flux band

spec_emin = 0.1  # keV, minimum energy of mock spectrum
spec_emax = 10.0  # keV, max energy of mock spectrum


def get_flux_scale(ind, fb_emin, fb_emax, spec_emin, spec_emax):
    f_g = np.log(spec_emax / spec_emin) * np.ones(ind.size)
    f_E = np.log(fb_emax / fb_emin) * np.ones(ind.size)
    n1 = ind != 1.0
    n2 = ind != 2.0
    f_g[n1] = (spec_emax ** (1.0 - ind[n1]) - spec_emin ** (1.0 - ind[n1])) / (
        1.0 - ind[n1]
    )
    f_E[n2] = (fb_emax ** (2.0 - ind[n2]) - fb_emin ** (2.0 - ind[n2])) / (
        2.0 - ind[n2]
    )
    fscale = f_g / f_E
    return fscale


def generate_fluxes(fov, prng):
    from soxs.data import cdf_agn, cdf_fluxes, cdf_gal

    prng = parse_prng(prng)

    fov = parse_value(fov, "arcmin")

    logf = np.log10(cdf_fluxes)

    n_gal = np.rint(cdf_gal[-1])
    n_agn = np.rint(cdf_agn[-1])
    F_gal = cdf_gal / cdf_gal[-1]
    F_agn = cdf_agn / cdf_agn[-1]
    f_gal = InterpolatedUnivariateSpline(F_gal, logf)
    f_agn = InterpolatedUnivariateSpline(F_agn, logf)

    fov_area = fov**2

    n_gal = int(n_gal * fov_area / 3600.0)
    n_agn = int(n_agn * fov_area / 3600.0)
    mylog.debug("%s AGN, %s galaxies in the FOV.", n_agn, n_gal)

    randvec1 = prng.uniform(size=n_agn)
    agn_fluxes = 10 ** f_agn(randvec1)

    randvec2 = prng.uniform(size=n_gal)
    gal_fluxes = 10 ** f_gal(randvec2)

    return agn_fluxes, gal_fluxes


def generate_positions(num, fov, sky_center, prng):
    fov *= 60.0  # convert to arcsec
    x = prng.uniform(low=-0.5 * fov, high=0.5 * fov, size=num)
    y = prng.uniform(low=-0.5 * fov, high=0.5 * fov, size=num)
    w = construct_wcs(*sky_center)
    ra0, dec0 = w.wcs_pix2world(x, y, 1)
    return ra0, dec0


def generate_sources(fov, sky_center, prng=None):
    r"""
    Make a catalog of point sources.

    Parameters
    ----------
    fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The field of view in arcminutes.
    sky_center : array-like
        The center RA, Dec of the field of view in degrees.
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    """
    prng = parse_prng(prng)

    fov = parse_value(fov, "arcmin")

    agn_fluxes, gal_fluxes = generate_fluxes(fov, prng)

    fluxes = np.concatenate([agn_fluxes, gal_fluxes])

    ind = np.concatenate(
        [get_agn_index(np.log10(agn_fluxes)), gal_index * np.ones(gal_fluxes.size)]
    )

    ra0, dec0 = generate_positions(fluxes.size, fov, sky_center, prng)

    return ra0, dec0, fluxes, ind


def make_ptsrc_background(
    exp_time,
    fov,
    sky_center,
    absorb_model=None,
    nH=None,
    area=40000.0,
    input_sources=None,
    output_sources=None,
    dump_fluxes_band=None,
    diffuse_unresolved=True,
    prng=None,
):
    r"""
    Make a point-source background.

    Parameters
    ----------
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time of the observation in seconds.
    fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The field of view in arcminutes.
    sky_center : array-like
        The center RA, Dec of the field of view in degrees.
    absorb_model : string, optional
        The absorption model to use, "wabs" or "tbabs".
        Defaults to the value in the SOXS configuration file.
    nH : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The hydrogen column in units of 10**22 atoms/cm**2.
        Defaults to the value in the SOXS configuration file.
    area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The effective area in cm**2. It must be large enough
        so that a sufficiently large sample is drawn for the
        ARF. Default: 40000.
    input_sources : string, optional
        If set to a filename, input the source positions, fluxes,
        and spectral indices from an ASCII table instead of generating
        them. Default: None
    output_sources : string, optional
        If set to a filename, output the properties of the sources
        within the field of view to a file. Default: None
    diffuse_unresolved : boolean, optional
        Add a diffuse component across the entire field of view to represent
        the unresolved flux from sources at very small fluxes. Default: True
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    """
    prng = parse_prng(prng)

    exp_time = parse_value(exp_time, "s")
    fov = parse_value(fov, "arcmin")
    if nH is not None:
        nH = parse_value(nH, "1.0e22*cm**-2")
    area = parse_value(area, "cm**2")
    if input_sources is None:
        ra0, dec0, fluxes, ind = generate_sources(fov, sky_center, prng=prng)
    else:
        mylog.info("Reading in point-source properties from %s.", input_sources)
        t = ascii.read(input_sources)
        ra0 = t["RA"].data
        dec0 = t["Dec"].data
        fluxes = t["flux_0.5_2.0_keV"].data
        ind = t["index"].data
    num_sources = fluxes.size

    mylog.debug("Generating spectra from %d sources.", num_sources)

    # If requested, output the source properties to a file
    if output_sources is not None:
        t = Table(
            [ra0, dec0, fluxes, ind], names=("RA", "Dec", "flux_0.5_2.0_keV", "index")
        )
        t["RA"].unit = "deg"
        t["Dec"].unit = "deg"
        t["flux_0.5_2.0_keV"].unit = "erg/(cm**2*s)"
        t["index"].unit = ""
        t.write(output_sources, format="ascii.ecsv", overwrite=True)

    # Add in the diffuse component for completely unresolved sources
    if diffuse_unresolved:
        """
        We add a diffuse background across the field of view, from completely
        unresolved sources, to match the fluxes in Table 3 of Hickox & Markevitch
        2006, ApJ, 645, 95. The first number "F12" is in erg/s/cm**2/deg**2, and
        is the difference between the total flux in sources below 5e-17 erg/s/cm**2
        and the "average" value in Table 3 for the unresolved flux in the 1-2 keV
        band of 1.04e-12 erg/s/cm**2/deg**2. The factor of 2.0 comes from the
        scaling between the fluxes in the 1-2 keV and 0.5-2 keV bands, assuming
        a spectral index of 2. The whole thing finally needs to be scaled by the
        FOV at the end.
        """
        F12 = 0.676e-12  # erg/s/cm**2/deg**2 in 1-2 keV band
        diffuse_flux = 2.0 * F12 * (fov / 60) ** 2
        fluxes = np.append(fluxes, diffuse_flux)
        ind = np.append(ind, 2.0)

    # Pre-calculate for optimization
    eratio = spec_emax / spec_emin
    oma = 1.0 - ind
    invoma = 1.0 / oma
    invoma[oma == 0.0] = 1.0
    fac1 = spec_emin**oma
    fac2 = spec_emax**oma - fac1
    fluxscale = get_flux_scale(ind, fb_emin, fb_emax, spec_emin, spec_emax)

    # Using the energy flux, determine the photon flux by simple scaling
    ref_ph_flux = fluxes * fluxscale * keV_per_erg
    # Now determine the number of photons we will generate
    n_photons = prng.poisson(ref_ph_flux * exp_time * area)

    all_energies = []
    all_ra = []
    all_dec = []
    detected = []

    for i, nph in enumerate(n_photons):

        if nph > 0:

            # Generate the energies in the source frame
            u = prng.uniform(size=nph)
            if ind[i] == 1.0:
                energies = spec_emin * (eratio**u)
            else:
                energies = fac1[i] + u * fac2[i]
                if invoma[i] == -1.0:
                    energies = 1.0 / energies
                else:
                    energies **= invoma[i]

            # Assign positions for this source
            if i < num_sources:
                ra = ra0[i] * np.ones(nph)
                dec = dec0[i] * np.ones(nph)
            else:
                # this is the diffuse source--fill the whole FOV
                ra, dec = generate_positions(nph, fov, sky_center, prng)

            all_energies.append(energies)
            all_ra.append(ra)
            all_dec.append(dec)
            detected.append(i)

    mylog.debug("Finished generating spectra.")

    if dump_fluxes_band is not None:
        dfluxes = np.zeros_like(ref_ph_flux)
        for idx, e in zip(detected, all_energies):
            eidxs = (e > dump_fluxes_band[0]) & (e < dump_fluxes_band[1])
            dfluxes[idx] = np.sum(e[eidxs]) * erg_per_keV / area / exp_time
        filename = f"dump_fluxes_{dump_fluxes_band[0]}_{dump_fluxes_band[1]}.dat"
        np.savetxt(filename, dfluxes)

    all_energies = np.concatenate(all_energies)
    all_ra = np.concatenate(all_ra)
    all_dec = np.concatenate(all_dec)

    # Remove some photons due to Galactic foreground absorption.
    # We will throw a lot of stuff away, but this is more general and still
    # faster.
    if nH is None:
        nH = float(soxs_cfg.get("soxs", "bkgnd_nH"))
    if nH > 0.0:
        if absorb_model is None:
            absorb_model = soxs_cfg.get("soxs", "bkgnd_absorb_model")
        if absorb_model == "wabs":
            absorb = get_wabs_absorb(all_energies, nH)
        elif absorb_model == "tbabs":
            absorb = get_tbabs_absorb(all_energies, nH)
        randvec = prng.uniform(size=all_energies.size)
        all_energies = all_energies[randvec < absorb]
        all_ra = all_ra[randvec < absorb]
        all_dec = all_dec[randvec < absorb]
        all_nph = all_energies.size
        mylog.debug("%d photons remain after foreground galactic absorption.", all_nph)

    all_flux = np.sum(all_energies) * erg_per_keV / (exp_time * area)

    output_events = {
        "ra": all_ra,
        "dec": all_dec,
        "energy": all_energies,
        "flux": all_flux,
    }

    return output_events


def make_point_sources_file(
    filename,
    name,
    exp_time,
    fov,
    sky_center,
    absorb_model=None,
    nH=None,
    area=40000.0,
    append=False,
    overwrite=False,
    src_filename=None,
    input_sources=None,
    output_sources=None,
    diffuse_unresolved=True,
    prng=None,
):
    """
    Make a SIMPUT catalog made up of contributions from
    point sources.

    Parameters
    ----------
    filename : string
        The filename for the SIMPUT catalog.
    name : string
        The name of the SIMPUT photon list.
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time of the observation in seconds.
    fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The field of view in arcminutes.
    sky_center : array-like
        The center RA, Dec of the field of view in degrees.
    absorb_model : string, optional
        The absorption model to use, "wabs" or "tbabs".
        Defaults to the value in the SOXS configuration file.
    nH : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The hydrogen column in units of 10**22 atoms/cm**2.
        Defaults to the value in the SOXS configuration file.
    area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The effective area in cm**2. It must be large enough
        so that a sufficiently large sample is drawn for the
        ARF. Default: 40000.
    append : boolean, optional
        If True, the photon list source will be appended to an existing
        SIMPUT catalog. Default: False
    overwrite : boolean, optional
        Set to True to overwrite previous files. Default: False
    src_filename : string, optional
        If set, this will be the filename to write the source
        to. By default, the source will be written to the same
        file as the SIMPUT catalog.
    input_sources : string, optional
        If set to a filename, input the source positions, fluxes,
        and spectral indices from an ASCII table instead of generating
        them. Default: None
    output_sources : string, optional
        If set to a filename, output the properties of the sources
        within the field of view to a file. Default: None
    diffuse_unresolved : boolean, optional
        Add a diffuse component across the entire field of view to represent
        the unresolved flux from sources at very small fluxes. Default: True
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    """
    events = make_ptsrc_background(
        exp_time,
        fov,
        sky_center,
        absorb_model=absorb_model,
        nH=nH,
        area=area,
        input_sources=input_sources,
        output_sources=output_sources,
        prng=prng,
        diffuse_unresolved=diffuse_unresolved,
    )
    phlist = SimputPhotonList(
        events["ra"], events["dec"], events["energy"], events["flux"], name=name
    )
    if append:
        cat = SimputCatalog.from_file(filename)
        cat.append(phlist, src_filename=src_filename, overwrite=overwrite)
    else:
        cat = SimputCatalog.from_source(
            filename, phlist, src_filename=src_filename, overwrite=overwrite
        )
    return cat


def make_point_source_list(output_file, fov, sky_center, prng=None):
    r"""
    Make a list of point source properties and write it to an ASCII
    table file.

    Parameters
    ----------
    output_file : string
        The ASCII table file to write the source properties to.
    fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The field of view in arcminutes.
    sky_center : array-like
        The center RA, Dec of the field of view in degrees.
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    """
    ra0, dec0, fluxes, ind = generate_sources(fov, sky_center, prng=prng)

    t = Table(
        [ra0, dec0, fluxes, ind], names=("RA", "Dec", "flux_0.5_2.0_keV", "index")
    )
    t["RA"].unit = "deg"
    t["Dec"].unit = "deg"
    t["flux_0.5_2.0_keV"].unit = "erg/(cm**2*s)"
    t["index"].unit = ""
    t.write(output_file, format="ascii.ecsv", overwrite=True)
