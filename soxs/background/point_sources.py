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

    cdf_fluxes = np.insert(np.logspace(-19, -13, 35)[::-1], 0, 1.0)

    cdf_ngal = np.array([0.0, 1.10390e-01, 1.68099e-01, 2.70088e-01,
                         4.19863e-01, 6.41722e-01, 1.00641e+00,
                         1.64646e+00, 2.56843e+00, 4.03708e+00,
                         6.44180e+00, 1.04667e+01, 1.67192e+01,
                         2.80757e+01, 4.58018e+01, 7.47671e+01,
                         1.27418e+02, 2.16680e+02, 3.57036e+02,
                         6.11174e+02, 1.08117e+03, 1.90768e+03,
                         3.28902e+03, 5.71005e+03, 9.65142e+03,
                         1.55802e+04, 2.40789e+04, 3.61239e+04,
                         5.16013e+04, 7.20866e+04, 9.68431e+04,
                         1.27856e+05, 1.64685e+05, 2.09966e+05,
                         2.63367e+05, 3.28799e+05])

    bl_nagn = np.array([0.0, 2.358680e+00, 5.798860e+00, 9.015150e+00,
                        2.206960e+01, 3.399270e+01, 7.703720e+01,
                        1.152240e+02, 2.200780e+02, 3.095730e+02,
                        5.098220e+02, 6.748330e+02, 1.005770e+03,
                        1.271110e+03, 1.767730e+03, 2.156880e+03,
                        2.876720e+03, 3.431040e+03, 4.460510e+03,
                        5.240910e+03, 6.610130e+03, 7.628250e+03,
                        9.299180e+03, 1.051860e+04, 1.248300e+04,
                        1.389550e+04, 1.605170e+04, 1.757760e+04,
                        1.982030e+04, 2.138530e+04, 2.373840e+04,
                        2.536190e+04, 2.776310e+04, 2.939940e+04,
                        3.048940e+04, 3.121540e+04])

    av_nagn = np.array([0.0, 2.597383e+00, 4.876015e+00, 9.153653e+00,
                        1.718399e+01, 3.225919e+01, 5.724188e+01,
                        9.271571e+01, 1.430281e+02, 2.143854e+02,
                        3.155914e+02, 4.582297e+02, 6.388914e+02,
                        8.593574e+02, 1.128393e+03, 1.456700e+03,
                        1.857341e+03, 2.346243e+03, 2.942856e+03,
                        3.670910e+03, 4.559362e+03, 5.643552e+03,
                        6.966600e+03, 8.581144e+03, 1.055136e+04,
                        1.295565e+04, 1.588962e+04, 1.947002e+04,
                        2.383920e+04, 2.917090e+04, 3.567732e+04,
                        4.361714e+04, 5.330617e+04, 6.512993e+04,
                        7.955844e+04, 9.716574e+04])

    if cdf_type == "av":
        cdf_nagn = av_nagn
    else:
        cdf_nagn = bl_nagn

    n_gal = np.rint(cdf_ngal[-1])
    n_agn = np.rint(cdf_nagn[-1])
    cdf_ngal /= cdf_ngal[-1]
    cdf_nagn /= cdf_nagn[-1]

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

    n_gal = int(n_gal*fov_area/3600.0)
    n_agn = int(n_agn*fov_area/3600.0)
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
