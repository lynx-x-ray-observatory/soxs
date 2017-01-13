import os
import numpy as np

from astropy.cosmology import FlatLambdaCDM

from soxs.simput import write_photon_list
from soxs.spatial import BetaModel, construct_wcs
from soxs.spectra import ApecGenerator
from soxs.utils import soxs_files_path, mylog

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

lum_file = os.path.join(soxs_files_path, "lum_table.txt")
halos_cat_file = os.path.join(soxs_files_path, "halo_catalog.dat")

def lum(M_mean, z_mean): 
    # Alexey's cosmology II paper, eq. 22
    E_z = np.sqrt(omega_m * (1 + z_mean) ** 3 + omega_k * (1 + z_mean) ** 2 + omega_l)
    ln_Lx = 47.392 + 1.61*np.log(M_mean) + 1.850 * np.log(E_z) - 0.39*np.log(h0/0.72)
    return np.exp(ln_Lx)

def Tx(M_mean, z_mean): 
    # Alexey's cosmology II paper, eq. 6, Table 3 for coefficients.
    E_z = np.sqrt(omega_m * (1 + z_mean) ** 3 + omega_k * (1 + z_mean) ** 2 + omega_l)
    return 5.0 * (M_mean*E_z*h0/3.02e14)**(1.0/1.53)

def flux2lum(kT,z):
    if kT > 15.0: kT = 15.0 ## to prevent the program from crashing in case a very massive halo happens to be on the FOV
    line = round((kT-0.1)/0.1)
    column = round(z/0.05)
    flux2lum = lum_table[line][column]
    return flux2lum

def make_cosmological_background(simput_prefix, phlist_prefix, exp_time, fov, sky_center,
                                 nH=0.05, area=40000.0, append=False, clobber=False, 
                                 prng=np.random):

    cosmo = FlatLambdaCDM(H0=100.0*h0, Om0=omega_m)
    agen = ApecGenerator(0.1, 10.0, 10000)

    mylog.info("Loading halo data from catalog: %s" % halos_cat_file)
    halo_data = np.loadtxt(halos_cat_file)

    scale = cosmo.kpc_proper_per_arcmin(halo_data[:,0]).to("Mpc/arcmin")

    # Scaling masses and transverse distances by hubble parameter,
    # converting transverse distances to arcmin
    halo_data[:,1:] /= h0
    halo_data[:,2:] /= scale.value

    # 600. arcmin = 10 degrees (total FOV of catalog = 100 deg^2)
    fov_cat = 10.0*60.0
    w = construct_wcs(*sky_center)

    xc, yc = prng.uniform(low=0.5*(fov-fov_cat), high=0.5*(fov_cat-fov), size=2)
    xlo = xc-1.1*0.5*fov
    xhi = xc+1.1*0.5*fov
    ylo = yc-1.1*0.5*fov
    yhi = yc+1.1*0.5*fov

    mylog.info("Selecting halos in the FOV.")

    fov_idxs = (halo_data[:,2] >= xlo) & (halo_data[:,2] <= xhi)
    fov_idxs = (halo_data[:,3] >= ylo) & (halo_data[:,3] <= yhi) & fov_idxs

    n_halos = fov_idxs.sum()

    mylog.info("Number of halos in the field of view: %d" % n_halos)

    # Now select the specific halos which are in the FOV
    z = halo_data[fov_idxs,0]
    m = halo_data[fov_idxs,1]
    s = scale[fov_idxs].to("kpc/arcsec").value
    ra0, dec0 = w.wcs_pix2world(halo_data[fov_idxs,2], halo_data[fov_idxs,3], 1)

    # Some cosmological stuff
    rho_crit = cosmo.critical_density(z).to("Msun/kpc**3").value

    # halo temperature and k-corrected flux
    kT = Tx(m, z)
    flux_kcorr = 1.0e-14*lum(m, z)/flux2lum(kT, z)

    # halo scale radius
    r500 = (3.0*m/(4.0*np.pi*500*rho_crit))**(1.0/3.0)
    rc = r500/conc/s

    # Halo slope parameter
    beta = prng.normal(loc=0.666, scale=0.05, size=n_halos)
    beta[beta < 0.5] = 0.5

    # Halo ellipticity
    ellip = prng.normal(loc=0.85, scale=0.15, size=n_halos)
    ellip[ellip < 0.0] = 1.0e-3

    # Halo orientation
    theta = 2.0*np.pi*prng.uniform(size=n_halos)

    tot_flux = 0.0
    ee = []
    ra = []
    dec = []

    for halo in range(n_halos):
        spec = agen.get_spectrum(kT[halo], abund, z[halo], 1.0)
        spec.rescale_flux(flux_kcorr[halo], emin=emin, emax=emax, flux_type="energy")
        spec.apply_foreground_absorption(nH)
        e = spec.generate_energies(exp_time, area, prng=prng)
        beta_model = BetaModel(ra0[halo], dec0[halo], rc[halo], beta[halo], e.size, 
                               ellipticity=ellip[halo], theta=theta[halo])
        tot_flux += e.flux
        ee.append(e.value)
        ra.append(beta_model.ra.value)
        ra.append(beta_model.dec.value)

    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    ee = np.concatenate(ee)

    write_photon_list(simput_prefix, phlist_prefix, tot_flux, ra, dec, ee, 
                      append=append, clobber=clobber)

