import os

import h5py
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from tqdm.auto import tqdm

from soxs.simput import SimputCatalog, SimputPhotonList
from soxs.spatial import BetaModel, construct_wcs
from soxs.thermal_spectra import ApecGenerator
from soxs.utils import mylog, parse_prng, parse_value, soxs_files_path

# Cosmological parameters for the catalog
# SHOULD NOT BE ALTERED
omega_m = 0.279
omega_l = 0.721
omega_k = 1.0 - omega_m - omega_l
h0 = 0.7

# Parameters for the halo fluxes
# SHOULD NOT BE ALTERED
emin = 0.5
emax = 2.0
abund = 0.3
conc = 10.0

lum_table_file = os.path.join(soxs_files_path, "lum_table.h5")
halos_cat_file = os.path.join(soxs_files_path, "halo_catalog.h5")


def lum(M_mean, z_mean):
    # Alexey's cosmology II paper, eq. 22
    E_z = np.sqrt(
        omega_m * (1.0 + z_mean) ** 3 + omega_k * (1.0 + z_mean) ** 2 + omega_l
    )
    ln_Lx = (
        47.392 + 1.61 * np.log(M_mean) + 1.850 * np.log(E_z) - 0.39 * np.log(h0 / 0.72)
    )
    return np.exp(ln_Lx)


def Tx(M_mean, z_mean):
    # Alexey's cosmology II paper, eq. 6, Table 3 for coefficients.
    E_z = np.sqrt(
        omega_m * (1.0 + z_mean) ** 3 + omega_k * (1.0 + z_mean) ** 2 + omega_l
    )
    return 5.0 * (M_mean * E_z * h0 / 3.02e14) ** (1.0 / 1.53)


def flux2lum(kT, z):
    lum_table = h5py.File(lum_table_file, "r")
    kT_idxs = np.round((kT - 0.1) / 0.1).astype("int")
    z_idxs = np.round(z / 0.05).astype("int")
    flux2lum = lum_table["Lx"][()][kT_idxs, z_idxs]
    lum_table.close()
    return flux2lum


def make_cosmological_sources(
    exp_time,
    fov,
    sky_center,
    cat_center=None,
    absorb_model="wabs",
    nH=0.05,
    area=40000.0,
    output_sources=None,
    write_regions=None,
    prng=None,
):
    r"""
    Make an X-ray source made up of contributions from
    galaxy clusters, galaxy groups, and galaxies.

    Parameters
    ----------
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time of the observation in seconds.
    fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The field of view in arcminutes.
    sky_center : array-like
        The center RA, Dec of the field of view in degrees.
    cat_center : array-like
        The center of the field in the coordinates of the
        halo catalog, which range from -5.0 to 5.0 in
        degrees in both directions. If None is given, a
        center will be randomly chosen.
    absorb_model : string, optional
        The absorption model to use, "wabs" or "tbabs". Default: "wabs"
    nH : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The hydrogen column in units of 10**22 atoms/cm**2.
        Default: 0.05
    area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The effective area in cm**2. It must be large enough
        so that a sufficiently large sample is drawn for the
        ARF. Default: 40000.
    output_sources : string, optional
        If set to a filename, output the properties of the sources
        within the field of view to a file. Default: None
    write_regions : string, optional
        If set to a filename, output circle ds9 regions corresponding to the
        positions of the halos with radii corresponding to their R500
        projected on the sky.  Default: None
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    """
    exp_time = parse_value(exp_time, "s")
    fov = parse_value(fov, "arcmin")
    if nH is not None:
        nH = parse_value(nH, "1.0e22*cm**-2")
    area = parse_value(area, "cm**2")
    prng = parse_prng(prng)
    cosmo = FlatLambdaCDM(H0=100.0 * h0, Om0=omega_m)
    agen = ApecGenerator(0.1, 10.0, 10000, broadening=False)

    mylog.info("Creating photons from cosmological sources.")

    mylog.info("Loading halo data from catalog: %s", halos_cat_file)
    halo_data = h5py.File(halos_cat_file, "r")

    scale = cosmo.kpc_comoving_per_arcmin(halo_data["redshift"][()]).to("Mpc/arcmin")

    # 600. arcmin = 10 degrees (total FOV of catalog = 100 deg^2)
    fov_cat = 10.0 * 60.0
    w = construct_wcs(*sky_center)

    cat_min = -0.5 * fov_cat
    cat_max = 0.5 * fov_cat

    if cat_center is None:
        xc, yc = prng.uniform(low=cat_min + 0.5 * fov, high=cat_max - 0.5 * fov, size=2)
    else:
        xc, yc = cat_center
        xc *= 60.0
        yc *= 60.0
        xc, yc = np.clip([xc, yc], cat_min + 0.5 * fov, cat_max - 0.5 * fov)

    mylog.info(
        "Coordinates of the FOV within the catalog are (%g, %g) deg.",
        xc / 60.0,
        yc / 60.0,
    )

    xlo = xc - 1.1 * 0.5 * fov
    xhi = xc + 1.1 * 0.5 * fov
    ylo = yc - 1.1 * 0.5 * fov
    yhi = yc + 1.1 * 0.5 * fov

    mylog.info("Selecting halos in the FOV.")

    halo_x = halo_data["x"][()].astype("float64") / (h0 * scale.value)
    halo_y = halo_data["y"][()].astype("float64") / (h0 * scale.value)

    fov_idxs = (halo_x >= xlo) & (halo_x <= xhi)
    fov_idxs = (halo_y >= ylo) & (halo_y <= yhi) & fov_idxs

    n_halos = fov_idxs.sum()

    mylog.info("Number of halos in the field of view: %d", n_halos)

    # Now select the specific halos which are in the FOV
    z = halo_data["redshift"][fov_idxs].astype("float64")
    m = halo_data["M500c"][fov_idxs].astype("float64") / h0
    # We need to compute proper scales here
    s = scale[fov_idxs].to("Mpc/arcsec").value / (1.0 + z)
    ra0, dec0 = w.wcs_pix2world(
        (halo_x[fov_idxs] - xc) * 60.0, (halo_y[fov_idxs] - yc) * 60.0, 1
    )

    # Close the halo catalog file
    halo_data.close()

    # Some cosmological stuff
    rho_crit = cosmo.critical_density(z).to("Msun/Mpc**3").value

    # halo temperature and k-corrected flux
    kT = Tx(m, z)
    flux_kcorr = 1.0e-14 * lum(m, z) / flux2lum(kT, z)

    # halo scale radius
    r500 = (3.0 * m / (4.0 * np.pi * 500 * rho_crit)) ** (1.0 / 3.0)
    r500_kpc = r500 * 1000.0
    rc_kpc = r500 / conc * 1000.0
    rc = r500 / conc / s

    # Halo slope parameter
    beta = prng.normal(loc=0.666, scale=0.05, size=n_halos)
    beta[beta < 0.5] = 0.5

    # Halo ellipticity
    ellip = prng.normal(loc=0.85, scale=0.15, size=n_halos)
    ellip[ellip < 0.0] = 1.0e-3

    # Halo orientation
    theta = 360.0 * prng.uniform(size=n_halos)

    # If requested, output the source properties to a file
    if output_sources is not None:
        t = Table(
            [ra0, dec0, rc_kpc, beta, ellip, theta, m, r500_kpc, kT, z, flux_kcorr],
            names=(
                "RA",
                "Dec",
                "r_c",
                "beta",
                "ellipticity",
                "theta",
                "M500c",
                "r500",
                "kT",
                "redshift",
                "flux_0.5_2.0_keV",
            ),
        )
        t["RA"].unit = "deg"
        t["Dec"].unit = "deg"
        t["flux_0.5_2.0_keV"].unit = "erg/(cm**2*s)"
        t["r_c"].unit = "kpc"
        t["theta"].unit = "deg"
        t["M500c"].unit = "solMass"
        t["r500"].unit = "kpc"
        t["kT"].unit = "kT"
        t.write(output_sources, format="ascii.ecsv", overwrite=True)

    if write_regions is not None:
        from astropy.coordinates import Angle, SkyCoord
        from regions import CircleSkyRegion, write_ds9

        regs = []
        for halo in range(n_halos):
            c = SkyCoord(ra0[halo], dec0[halo], unit=("deg", "deg"), frame="fk5")
            scale = cosmo.kpc_proper_per_arcmin(z[halo]).to("kpc/deg")
            r500c = r500_kpc / scale.value
            r = Angle(r500c, "deg")
            reg = CircleSkyRegion(c, r)
            regs.append(reg)
        write_ds9(regs, write_regions)

    tot_flux = 0.0
    ee = []
    ra = []
    dec = []

    pbar = tqdm(leave=True, total=n_halos, desc="Generating photons from halos ")
    for halo in range(n_halos):
        spec = agen.get_spectrum(kT[halo], abund, z[halo], 1.0)
        spec.rescale_flux(flux_kcorr[halo], emin=emin, emax=emax, flux_type="energy")
        if nH is not None:
            spec.apply_foreground_absorption(nH, model=absorb_model)
        e = spec.generate_energies(exp_time, area, prng=prng, quiet=True)
        beta_model = BetaModel(
            ra0[halo],
            dec0[halo],
            rc[halo],
            beta[halo],
            ellipticity=ellip[halo],
            theta=theta[halo],
        )
        xsky, ysky = beta_model.generate_coords(e.size, prng=prng)
        tot_flux += e.flux
        ee.append(e.value)
        ra.append(xsky.value)
        dec.append(ysky.value)
        pbar.update()
    pbar.close()

    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    ee = np.concatenate(ee)

    mylog.info("Created %d photons from cosmological sources.", ee.size)

    output_events = {"ra": ra, "dec": dec, "energy": ee, "flux": tot_flux.value}

    return output_events


def make_cosmological_sources_file(
    filename,
    name,
    exp_time,
    fov,
    sky_center,
    cat_center=None,
    absorb_model="wabs",
    nH=0.05,
    area=40000.0,
    overwrite=False,
    output_sources=None,
    write_regions=None,
    src_filename=None,
    prng=None,
    append=False,
):
    r"""
    Make a SIMPUT catalog made up of contributions from
    galaxy clusters, galaxy groups, and galaxies.

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
    cat_center : array-like
        The center of the field in the coordinates of the
        halo catalog, which range from -5.0 to 5.0 degrees
        along both axes. If None is given, a center will be
        randomly chosen.
    absorb_model : string, optional
        The absorption model to use, "wabs" or "tbabs". Default: "wabs"
    nH : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The hydrogen column in units of 10**22 atoms/cm**2.
        Default: 0.05
    area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The effective area in cm**2. It must be large enough
        so that a sufficiently large sample is drawn for the
        ARF. Default: 40000.
    overwrite : boolean, optional
        Set to True to overwrite previous files. Default: False
    output_sources : string, optional
        If set to a filename, output the properties of the sources
        within the field of view to an ASCII file. Default: None
    write_regions : string, optional
        If set to a filename, output circle ds9 regions corresponding to the
        positions of the halos with radii corresponding to their R500
        projected on the sky. Default: None
    src_filename : string, optional
        If set, this will be the filename to write the source
        to. By default, the source will be written to the same
        file as the SIMPUT catalog
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    append : boolean, optional
        If True, the photon list source will be appended to an existing
        SIMPUT catalog. Default: False
    """
    events = make_cosmological_sources(
        exp_time,
        fov,
        sky_center,
        cat_center=cat_center,
        absorb_model=absorb_model,
        nH=nH,
        area=area,
        output_sources=output_sources,
        write_regions=write_regions,
        prng=prng,
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
