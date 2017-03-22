import os
import numpy as np
import h5py

from tqdm import tqdm

from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

from soxs.simput import write_photon_list
from soxs.spatial import BetaModel, construct_wcs
from soxs.spectra import ApecGenerator
from soxs.utils import soxs_files_path, mylog, parse_prng

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
    E_z = np.sqrt(omega_m * (1 + z_mean) ** 3 + omega_k * (1 + z_mean) ** 2 + omega_l)
    ln_Lx = 47.392 + 1.61*np.log(M_mean) + 1.850 * np.log(E_z) - 0.39*np.log(h0/0.72)
    return np.exp(ln_Lx)

def Tx(M_mean, z_mean): 
    # Alexey's cosmology II paper, eq. 6, Table 3 for coefficients.
    E_z = np.sqrt(omega_m * (1 + z_mean) ** 3 + omega_k * (1 + z_mean) ** 2 + omega_l)
    return 5.0 * (M_mean*E_z*h0/3.02e14)**(1.0/1.53)

def flux2lum(kT, z):
    lum_table = h5py.File(lum_table_file, "r")
    kT_idxs = np.round((kT-0.1)/0.1).astype('int')
    z_idxs = np.round(z/0.05).astype('int')
    flux2lum = lum_table["Lx"].value[kT_idxs, z_idxs]
    lum_table.close()
    return flux2lum

def make_cosmological_sources(exp_time, fov, sky_center, cat_center=None,
                              nH=0.05, area=40000.0, output_sources=None,
                              prng=None):
    r"""
    Make an X-ray source made up of contributions from
    galaxy clusters, galaxy groups, and galaxies. 

    Parameters
    ----------
    exp_time : float
        The exposure time of the observation in seconds.
    fov : float
        The field of view in arcminutes.
    sky_center : array-like
        The center RA, Dec of the field of view in degrees.
    cat_center : array-like
        The center of the field in the coordinates of the
        halo catalog, which range from -5.0 to 5.0 in 
        degrees in both directions. If None is given, a 
        center will be randomly chosen.
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
    cosmo = FlatLambdaCDM(H0=100.0*h0, Om0=omega_m)
    agen = ApecGenerator(0.1, 10.0, 10000, broadening=False)

    mylog.info("Creating photons from cosmological sources.")

    mylog.info("Loading halo data from catalog: %s" % halos_cat_file)
    halo_data = h5py.File(halos_cat_file, "r")

    scale = cosmo.kpc_proper_per_arcmin(halo_data["redshift"]).to("Mpc/arcmin")

    # 600. arcmin = 10 degrees (total FOV of catalog = 100 deg^2)
    fov_cat = 10.0*60.0
    w = construct_wcs(*sky_center)

    cat_min = -0.5*fov_cat
    cat_max = 0.5*fov_cat

    if cat_center is None:
        xc, yc = prng.uniform(low=cat_min+0.5*fov, high=cat_max-0.5*fov, size=2)
    else:
        xc, yc = cat_center
        xc *= 60.0
        yc *= 60.0
        xc, yc = np.clip([xc, yc], cat_min+0.5*fov, cat_max-0.5*fov)

    mylog.info("Coordinates of the FOV within the catalog are (%g, %g) deg." %
               (xc/60.0, yc/60.0))

    xlo = (xc-1.1*0.5*fov)*scale.value*h0
    xhi = (xc+1.1*0.5*fov)*scale.value*h0
    ylo = (yc-1.1*0.5*fov)*scale.value*h0
    yhi = (yc+1.1*0.5*fov)*scale.value*h0

    mylog.info("Selecting halos in the FOV.")

    fov_idxs = (halo_data["x"] >= xlo) & (halo_data["x"] <= xhi)
    fov_idxs = (halo_data["y"] >= ylo) & (halo_data["y"] <= yhi) & fov_idxs

    n_halos = fov_idxs.sum()

    mylog.info("Number of halos in the field of view: %d" % n_halos)

    # Now select the specific halos which are in the FOV
    z = halo_data["redshift"][fov_idxs].astype("float64")
    m = halo_data["M500c"][fov_idxs].astype("float64")/h0
    s = scale[fov_idxs].to("Mpc/arcsec").value
    ra0, dec0 = w.wcs_pix2world(halo_data["x"][fov_idxs]/(h0*s)-xc*60.0,
                                halo_data["y"][fov_idxs]/(h0*s)-yc*60.0, 1)

    # Close the halo catalog file
    halo_data.close()

    # Some cosmological stuff
    rho_crit = cosmo.critical_density(z).to("Msun/Mpc**3").value

    # halo temperature and k-corrected flux
    kT = Tx(m, z)
    flux_kcorr = 1.0e-14*lum(m, z)/flux2lum(kT, z)

    # halo scale radius
    r500 = (3.0*m/(4.0*np.pi*500*rho_crit))**(1.0/3.0)
    r500_kpc = r500 * 1000.0
    rc_kpc = r500/conc * 1000.0
    rc = r500/conc/s

    # Halo slope parameter
    beta = prng.normal(loc=0.666, scale=0.05, size=n_halos)
    beta[beta < 0.5] = 0.5

    # Halo ellipticity
    ellip = prng.normal(loc=0.85, scale=0.15, size=n_halos)
    ellip[ellip < 0.0] = 1.0e-3

    # Halo orientation
    theta = 360.0*prng.uniform(size=n_halos)

    # If requested, output the source properties to a file
    if output_sources is not None:
        t = Table([ra0, dec0, rc_kpc, beta, ellip, theta, m, r500_kpc, kT, z, flux_kcorr],
                  names=('RA', 'Dec', 'r_c', 'beta', 'ellipticity',
                         'theta', 'M500c', 'r500', 'kT', 'redshift', 'flux_0.5_2.0_keV'))
        t["RA"].unit = "deg"
        t["Dec"].unit = "deg"
        t["flux_0.5_2.0_keV"].unit = "erg/(cm**2*s)"
        t["r_c"].unit = "kpc"
        t["theta"].unit = "deg"
        t["M500c"].unit = "solMass"
        t["r500"].unit = "kpc"
        t["kT"].unit = "kT"
        t.write(output_sources, format='ascii.ecsv', overwrite=True)

    tot_flux = 0.0
    ee = []
    ra = []
    dec = []

    pbar = tqdm(leave=True, total=n_halos, desc="Generating photons from halos ")
    for halo in range(n_halos):
        spec = agen.get_spectrum(kT[halo], abund, z[halo], 1.0)
        spec.rescale_flux(flux_kcorr[halo], emin=emin, emax=emax, flux_type="energy")
        if nH is not None:
            spec.apply_foreground_absorption(nH)
        e = spec.generate_energies(exp_time, area, prng=prng, quiet=True)
        beta_model = BetaModel(ra0[halo], dec0[halo], rc[halo], beta[halo], e.size, 
                               ellipticity=ellip[halo], theta=theta[halo], prng=prng)
        tot_flux += e.flux
        ee.append(e.value)
        ra.append(beta_model.ra.value)
        dec.append(beta_model.dec.value)
        pbar.update()
    pbar.close()

    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    ee = np.concatenate(ee)

    mylog.info("Created %d photons from cosmological sources." % ee.size)

    output_events = {"ra": ra, "dec": dec, "energy": ee, 
                     "flux": tot_flux.value}

    return output_events

def make_cosmological_sources_file(simput_prefix, phlist_prefix, exp_time, fov, 
                                   sky_center, cat_center=None, nH=0.05, 
                                   area=40000.0, append=False, overwrite=False, 
                                   output_sources=None, prng=None):
    r"""
    Make a SIMPUT catalog made up of contributions from
    galaxy clusters, galaxy groups, and galaxies.

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
    cat_center : array-like
        The center of the field in the coordinates of the
        halo catalog, which range from -5.0 to 5.0 degrees
        along both axes. If None is given, a center will be
        randomly chosen.
    nH : float, optional
        The hydrogen column in units of 10**22 atoms/cm**2. 
        Default: 0.05
    area : float, optional
        The effective area in cm**2. It must be large enough 
        so that a sufficiently large sample is drawn for the 
        ARF. Default: 40000.
    append : boolean, optional
        If True, append a new source an existing SIMPUT 
        catalog. Default: False
    overwrite : boolean, optional
        Set to True to overwrite previous files. Default: False
    output_sources : string, optional
        If set to a filename, output the properties of the sources
        within the field of view to a file. Default: None
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    events = make_cosmological_sources(exp_time, fov, sky_center, cat_center=cat_center,
                                       nH=nH, area=area, output_sources=output_sources,
                                       prng=prng)
    write_photon_list(simput_prefix, phlist_prefix, events["flux"], events["ra"], 
                      events["dec"], events["energy"], append=append, overwrite=overwrite)
