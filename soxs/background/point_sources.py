"""
Draw point sources from a logN-logS distribution and generate 
simulated events from them. Flux units are 10^-14 erg/cm^2/s.  
Author: Scott Randall (srandall@cfa.harvard.edu)
"""
import numpy as np
from soxs import write_photon_list
from soxs.constants import keV_per_erg, erg_per_keV
from soxs.spectra import get_wabs_absorb
from soxs.utils import mylog

cdf_fluxes = np.array([1.00000e-19, 1.50131e-19, 2.25393e-19, 3.38386e-19,
                       5.08022e-19, 7.62699e-19, 1.14505e-18, 1.71907e-18,
                       2.58086e-18, 3.87468e-18, 5.81709e-18, 8.73326e-18,
                       1.31113e-17, 1.96842e-17, 2.95521e-17, 4.43669e-17,
                       6.66085e-17, 1.00000e-16, 1.50131e-16, 2.25393e-16,
                       3.38386e-16, 5.08022e-16, 7.62699e-16, 1.14505e-15,
                       1.71907e-15, 2.58086e-15, 3.87468e-15, 5.81709e-15,
                       8.73326e-15, 1.31113e-14, 1.96842e-14, 2.95521e-14,
                       4.43669e-14, 6.66085e-14, 1.00000e-13, 1.0])

cdf_ngal = np.array([328799.0, 263367.0, 209966.0, 164685.0, 127856.0,
                     96843.1, 72086.6, 51601.3, 36123.9, 24078.9,
                     15580.2, 9651.42, 5710.05, 3289.02, 1907.68,
                     1081.17, 611.174, 357.036, 216.680, 127.418,
                     74.7671, 45.8018, 28.0757, 16.7192, 10.4667,
                     6.44180, 4.03708, 2.56843, 1.64646, 1.00641,
                     0.641722, 0.419863, 0.270088, 0.168099, 0.110390])

bl_nagn = np.array([31215.4, 30489.4, 29399.4, 27763.1, 25361.9,
                    23738.4, 21385.3, 19820.3, 17577.6, 16051.7,
                    13895.5, 12483.0, 10518.6, 9299.18, 7628.25,
                    6610.13, 5240.91, 4460.51, 3431.04, 2876.72,
                    2156.88, 1767.73, 1271.11, 1005.77, 674.833,
                    509.822, 309.573, 220.078, 115.224, 77.0372,
                    33.9927, 22.0696, 9.01515, 5.79886, 2.35868])

av_nagn = np.array([97165.74, 79558.44, 65129.93, 53306.17, 43617.14,
                    35677.32, 29170.90, 23839.20, 19470.02, 15889.62,
                    12955.65, 10551.36, 8581.144, 6966.600, 5643.552,
                    4559.362, 3670.910, 2942.856, 2346.243, 1857.341,
                    1456.700, 1128.393, 859.3574, 638.8914, 458.2297,
                    315.5914, 214.3854, 143.0281, 92.71571, 57.24188,
                    32.25919, 17.18399, 9.153653, 4.876015, 2.597383])

class Bgsrc:
    def __init__(self, src_type, flux, z, ind):
        self.src_type = src_type
        self.flux = flux
        self.z = z
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

def make_ptsrc_background(simput_prefix, phlist_prefix, exp_time, fov, sky_center,
                          nH=0.05, nH_int=None, area=40000.0, append=False, 
                          clobber=False, cdf_type="av", prng=np.random):

    sources = []
    src_types = ['agn', 'gal']

    # parameters for making event file
    # for now, we're leaving these not user-configurable

    fb_emin = 0.5  # keV, low energy bound for the logN-logS flux band
    fb_emax = 2.0  # keV, high energy bound for the logN-logS flux band
    spec_emin = 0.1 # keV, minimum energy of mock spectrum
    spec_emax = 10.0 # keV, max energy of mock spectrum

    agn_ind = 1.2 # AGN photon index
    agn_z = 2.0 # AGN redshift
    gal_ind = 1.2 # galaxy photon index
    gal_z = 0.8 # galaxy redshift

    if cdf_type == "av":
        cdf_nagn = av_nagn
    else:
        cdf_nagn = bl_nagn

    n_gal = np.rint(cdf_ngal[0])
    n_agn = np.rint(cdf_nagn[0])
    cdf_fluxes = np.fliplr([cdf_fluxes])[0]
    cdf_ngal = np.fliplr([cdf_ngal])[0]
    cdf_nagn = np.fliplr([cdf_nagn])[0]
    cdf_ngal /= cdf_ngal[-1]
    cdf_nagn /= cdf_nagn[-1]
    cdf_ngal = np.insert(cdf_ngal, 0, 0.0)
    cdf_nagn = np.insert(cdf_nagn, 0, 0.0)

    redshifts = {"agn": agn_z, "gal": gal_z}
    indices = {"agn": agn_ind, "gal": gal_ind}
    fluxscale = {}
    for src_type in src_types:
        fluxscale[src_type] = get_flux_scale(indices[src_type], fb_emin, fb_emax,
                                             spec_emin, spec_emax)

    eph_mean_erg = 1.0*erg_per_keV

    S_min_obs = eph_mean_erg/(exp_time*area)
    mylog.debug("Flux of %g erg/cm^2/s gives roughly "
                "one photon during exposure." % S_min_obs)
    fov_area = fov**2

    n_gal = np.rint(n_gal*fov_area/3600.0)
    n_agn = np.rint(n_agn*fov_area/3600.0)
    num_sources = n_gal + n_agn
    mylog.info("%d AGN, %d galaxies in the FOV." % (n_agn, n_gal))

    randvec1 = prng.uniform(size=n_agn)
    agn_fluxes = np.interp(randvec1, cdf_nagn, cdf_fluxes)

    randvec2 = prng.uniform(size=n_gal)
    gal_fluxes = np.interp(randvec2, cdf_ngal, cdf_fluxes)

    for S in agn_fluxes:
        thissrc = Bgsrc("agn", S, redshifts["agn"], indices["agn"])
        sources.append(thissrc)
    for S in gal_fluxes:
        thissrc = Bgsrc("gal", S, redshifts["gal"], indices["gal"])
        sources.append(thissrc)

    mylog.info("Generating spectra from %d sources." % len(sources))
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

    u_src = prng.uniform(size=num_sources)

    for i, source in enumerate(sources):

        # Using the energy flux, determine the photon flux by simple scaling
        ref_ph_flux = source.flux*fluxscale[source.src_type]*keV_per_erg
        # Now determine the number of photons we will generate
        nph = np.modf(ref_ph_flux*exp_time*area)
        nph = np.int64(nph[1]) + np.int64(nph[0] >= u_src[i])

        if nph > 0:
            # Generate the energies in the source frame
            u = prng.uniform(size=nph)
            if source.ind == 1.0:
                energies = spec_emin*(eratio**u)
            else:
                energies = fac1[source.src_type] + u*fac2[source.src_type]
                energies **= invoma[source.src_type]
            # Here is where we apply intrinsic absorption for galaxies and agn.
            # Local galactic absorption is done at the end.
            if nH_int is not None and source.src_type in ["agn", "gal"]:
                absorb = get_wabs_absorb(energies*(1.0+source.z), nH_int)
                randvec = prng.uniform(size=energies.size)
                energies = energies[randvec < absorb]
            new_nph = energies.size
            # Assign positions for this source
            ra = prng.random()*np.ones(new_nph)*fov/(60.0*dec_scal) + ra_min
            dec = prng.random()*np.ones(new_nph)*fov/60.0 + dec_min

            all_energies.append(energies)
            all_ra.append(ra)
            all_dec.append(dec)

    mylog.info("Finished generating spectra.")

    all_energies = np.concatenate(all_energies)
    all_ra = np.concatenate(all_ra)
    all_dec = np.concatenate(all_dec)

    all_nph = all_energies.size

    mylog.info("Generated %d photons from point sources." % all_nph)

    # Remove some of the photons due to Galactic foreground absorption.
    # We will throw a lot of stuff away, but this is more general and still
    # faster. 
    absorb = get_wabs_absorb(all_energies, nH)
    randvec = prng.uniform(size=all_energies.size)
    all_energies = all_energies[randvec < absorb]
    all_ra = all_ra[randvec < absorb]
    all_dec = all_dec[randvec < absorb]
    all_nph = all_energies.size
    mylog.info("%d photons remain after foreground galactic absorption." % all_nph)

    all_flux = np.sum(all_energies)*erg_per_keV/(exp_time*area)

    write_photon_list(simput_prefix, phlist_prefix, all_flux, all_ra, all_dec, 
                      all_energies, clobber=clobber, append=append)
