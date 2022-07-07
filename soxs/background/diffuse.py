from soxs.instrument import perform_dither
from soxs.spectra import ConvolvedSpectrum
from soxs.thermal_spectra import ApecGenerator
from soxs.utils import parse_prng, \
    parse_value, mylog, create_region, \
    get_data_file, get_rot_mat, soxs_cfg
import numpy as np
from astropy.io import fits
from regions import PixCoord

"""
XSPEC model used to create the "default" foreground spectrum
  model  apec + wabs*apec
            0.099
                1
                0
          1.7e-06
            0.018
            0.225
                1
                0
          7.3e-07

XSPEC model used to create the "lem" foreground spectrum
  model  apec + tbabs*(bapec+bapec)
            0.099
                1
                0
          1.7e-06
            0.018
            0.225
                1
                0
              100
          7.3e-07
              0.7
                1
                0
              100
         8.76e-08 
"""


def make_frgnd_spectrum(arf, rmf):
    bkgnd_nH = float(soxs_cfg.get("soxs", "bkgnd_nH"))
    absorb_model = soxs_cfg.get("soxs", "bkgnd_absorb_model")
    frgnd_spec_model = soxs_cfg.get("soxs", "frgnd_spec_model")
    frgnd_velocity = float(soxs_cfg.get("soxs", "frgnd_velocity"))
    agen = ApecGenerator(rmf.ebins[0], rmf.ebins[-1], rmf.n_e,
                         broadening=True)
    spec = agen.get_spectrum(0.225, 1.0, 0.0, 7.3e-7, velocity=frgnd_velocity)
    if frgnd_spec_model == "halosat":
        spec += agen.get_spectrum(0.7, 1.0, 0.0, 8.76e-8, velocity=frgnd_velocity)
    spec.apply_foreground_absorption(bkgnd_nH, model=absorb_model)
    spec += agen.get_spectrum(0.099, 1.0, 0.0, 1.7e-6)
    spec.restrict_within_band(emin=0.1)
    spec = ConvolvedSpectrum.convolve(spec, arf)
    return spec


def read_instr_spectrum(filename, ext_area):
    """
    Read an instrumental background spectrum from
    a FITS PHA file.

    Parameters
    ----------
    filename : string
        The path to the file containing the spectrum.
    """
    fn = get_data_file(filename)
    with fits.open(fn) as f:
        hdu = f["SPECTRUM"]
        if "COUNTS" in hdu.data.names:
            count_rate = hdu.data["COUNTS"]/hdu.header["EXPOSURE"]
        elif "COUNT_RATE" in hdu.data.names:
            count_rate = hdu.data["COUNT_RATE"]
        elif "RATE" in hdu.data.names:
            count_rate = hdu.data["RATE"]
        else:
            raise RuntimeError("Cannot find a field for either "
                               "counts or count rate!")
        count_rate /= ext_area
    return count_rate


def generate_channel_spectrum(count_rate, t_exp, solid_angle, prng=None):
    """
    Generate photon energy channels from this diffuse
    background spectrum given an exposure time,
    effective area, and solid angle.

    Parameters
    ----------
    t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time in seconds.
    solid_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The solid angle in arcmin**2.
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    t_exp = parse_value(t_exp, "s")
    solid_angle = parse_value(solid_angle, "arcmin**2")
    prng = parse_prng(prng)
    fac = t_exp * solid_angle # Backgrounds are normalized to 1 arcmin**2
    return prng.poisson(lam=count_rate*fac).astype('int')


def make_diffuse_background(foreground, instr_bkgnd, inst_spec, event_params, 
                            arf, rmf, prng=None):
    from collections import defaultdict

    prng = parse_prng(prng)

    if foreground:
        mylog.info("Adding in astrophysical foreground.")
        fspec = rmf.convolve_spectrum(make_frgnd_spectrum(arf, rmf),
                                      event_params["exposure_time"],
                                      noisy=False, rate=True)

    if instr_bkgnd:
        mylog.info("Adding in instrumental background.")
        bkgnd_spec = inst_spec["bkgnd"]
        if isinstance(bkgnd_spec[0], str):
            nchips = len(event_params["chips"])
            bkgnd_spec = [bkgnd_spec]*nchips

    n_frgnd = 0
    n_instr = 0

    channels = (np.arange(rmf.n_ch) + rmf.cmin).astype("int32")

    bkg_events = defaultdict(list)
    pixel_area = (event_params["plate_scale"]*60.0)**2
    for i, chip in enumerate(event_params["chips"]):
        rtype = chip[0]
        args = chip[1:]
        r, bounds = create_region(rtype, args, 0.0, 0.0)
        sa = (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*pixel_area
        ncts = np.zeros(rmf.n_ch, dtype='int')
        if instr_bkgnd:
            ispec = read_instr_spectrum(bkgnd_spec[i][0], bkgnd_spec[i][1])
            icts = generate_channel_spectrum(ispec, event_params["exposure_time"], sa, prng=prng)
            ncts += icts
            n_instr += icts.sum()
        if foreground:
            fcts = generate_channel_spectrum(fspec, event_params["exposure_time"], sa, prng=prng)
            ncts += fcts
            n_frgnd += fcts.sum()
        chan = np.concatenate([channels[i] * np.ones(n)
                               for i, n in enumerate(ncts)])
        n_events = chan.size
        detx = prng.uniform(low=bounds[0], high=bounds[1], size=n_events)
        dety = prng.uniform(low=bounds[2], high=bounds[3], size=n_events)
        if rtype in ["Box", "Rectangle"]:
            thisc = slice(None, None, None)
            n_det = n_events
        else:
            thisc = r.contains(PixCoord(detx, dety))
            n_det = thisc.sum()
        ch = chan[thisc].astype('int')
        e = rmf.ch_to_eb(ch, prng=prng)
        bkg_events["energy"].append(e)
        bkg_events[rmf.chan_type].append(ch)
        bkg_events["detx"].append(detx[thisc])
        bkg_events["dety"].append(dety[thisc])
        bkg_events["chip_id"].append(i*np.ones(n_det))
    for key in bkg_events:
        bkg_events[key] = np.concatenate(bkg_events[key])

    if bkg_events["energy"].size == 0:
        raise RuntimeError("No background events were detected!!!")
    if foreground:
        mylog.info(f"Making {n_frgnd} events from the galactic foreground.")
    if instr_bkgnd:
        mylog.info(f"Making {n_instr} events from the instrumental background.")

    n_e = bkg_events["energy"].size

    bkg_events['time'] = prng.uniform(size=n_e, low=0.0,
                                      high=event_params["exposure_time"])

    x_offset, y_offset = perform_dither(bkg_events["time"],
                                        event_params["dither_params"])

    rot_mat = get_rot_mat(event_params["roll_angle"])

    det = np.array([bkg_events["detx"] + x_offset -
                    event_params["aimpt_coords"][0] -
                    event_params["aimpt_shift"][0],
                    bkg_events["dety"] + y_offset -
                    event_params["aimpt_coords"][1] -
                    event_params["aimpt_shift"][1]])
    pix = np.dot(rot_mat.T, det)

    bkg_events["xpix"] = pix[0, :] + event_params['pix_center'][0]
    bkg_events["ypix"] = pix[1, :] + event_params['pix_center'][1]

    return bkg_events
