import numpy as np
import os
import astropy.io.fits as pyfits
from soxs.utils import mylog, parse_prng

key_map = {"telescope": "TELESCOP",
           "mission": "MISSION",
           "instrument": "INSTRUME",
           "channel_type": "CHANTYPE",
           "nchan": "PHA_BINS"}

def add_background_from_file(events, event_params, bkg_file):
    f = pyfits.open(bkg_file)

    hdu = f["EVENTS"]

    sexp = event_params["exposure_time"]
    bexp = hdu.header["EXPOSURE"]

    if event_params["exposure_time"] > hdu.header["EXPOSURE"]:
        raise RuntimeError("The background file does not have sufficient exposure! Source "
                           "exposure time %g, background exposure time %g." % (sexp, bexp))

    for k1, k2 in key_map.items():
        if event_params[k1] != hdu.header[k2]:
            raise RuntimeError("'%s' keyword does not match! %s vs. %s" % (k1, event_params[k1],
                                                                           hdu.header[k2]))
    rmf1 = os.path.split(event_params["rmf"])[-1]
    rmf2 = hdu.header["RESPFILE"]
    arf1 = os.path.split(event_params["arf"])[-1]
    arf2 = hdu.header["ANCRFILE"]
    if rmf1 != rmf2:
        raise RuntimeError("RMFs do not match! %s vs. %s" % (rmf1, rmf2))
    if arf1 != arf2:
        raise RuntimeError("ARFs do not match! %s vs. %s" % (arf1, arf2))

    idxs = hdu.data["TIME"] < sexp

    mylog.info("Adding %d background events from %s." % (idxs.sum(), bkg_file))

    if event_params["roll_angle"] == hdu.header["ROLL_PNT"]:
        xpix = hdu.data["X"][idxs]
        ypix = hdu.data["Y"][idxs]
    else:
        roll_angle = np.deg2rad(event_params["roll_angle"])
        rot_mat = np.array([[np.sin(roll_angle), -np.cos(roll_angle)],
                            [-np.cos(roll_angle), -np.sin(roll_angle)]])
        xpix, ypix = np.dot(rot_mat, np.array([hdu.data["DETX"][idxs], 
                                               hdu.data["DETY"][idxs]]))
        xpix += hdu.header["TCRPX2"]
        ypix += hdu.header["TCRPX3"]

    all_events = {}
    for key in ["detx", "dety", "time", "ccd_id", event_params["channel_type"]]:
        all_events[key] = np.concatenate([events[key], hdu.data[key.upper()][idxs]])
    all_events["xpix"] = np.concatenate([events["xpix"], xpix])
    all_events["ypix"] = np.concatenate([events["ypix"], ypix])
    all_events["energy"] = np.concatenate([events["energy"], hdu.data["ENERGY"][idxs]/1000.0])

    f.close()

    return all_events

def make_uniform_background(energy, event_params, rmf, prng=None):
    from soxs.instrument import perform_dither
    import pyregion._region_filter as filter

    prng = parse_prng(prng)

    bkg_events = {}

    n_events = energy.size

    nx = event_params["num_pixels"]
    bkg_events["detx"] = prng.uniform(low=0.5, high=nx+0.5, size=n_events)
    bkg_events["dety"] = prng.uniform(low=0.5, high=nx+0.5, size=n_events)
    bkg_events["energy"] = energy

    if event_params["chips"] is None:
        bkg_events["chip_id"] = np.zeros(n_events, dtype='int')
    else:
        bkg_events["chip_id"] = -np.ones(n_events, dtype='int')
        for i, chip in enumerate(event_params["chips"]):
            thisc = np.ones(n_events, dtype='bool')
            for reg in chip["region"]:
                rtype = reg[0]
                mask = reg[1]
                args = reg[2:]
                r = getattr(filter, rtype)(*args)
                inside = r.inside(bkg_events["detx"], bkg_events["dety"])
                if not mask:
                    inside = np.logical_not(inside)
                thisc = np.logical_and(thisc, inside)
            bkg_events["chip_id"][thisc] = chip["id"]

    keep = bkg_events["chip_id"] > -1

    for key in bkg_events:
        bkg_events[key] = bkg_events[key][keep]

    n_e = bkg_events["energy"].size

    bkg_events['time'] = prng.uniform(size=n_e, low=0.0,
                                      high=event_params["exposure_time"])

    x_offset, y_offset = perform_dither(bkg_events["time"],
                                        event_params["dither_params"])

    bkg_events["detx"] -= event_params['det_center'][0]
    bkg_events["dety"] -= event_params['det_center'][1]

    roll_angle = np.deg2rad(event_params["roll_angle"])
    rot_mat = np.array([[np.sin(roll_angle), -np.cos(roll_angle)],
                        [-np.cos(roll_angle), -np.sin(roll_angle)]])

    det = np.array([bkg_events["detx"] + x_offset - event_params["aimpt_coords"][0],
                    bkg_events["dety"] + y_offset - event_params["aimpt_coords"][1]])
    pix = np.dot(rot_mat, det)

    bkg_events["xpix"] = pix[0, :] + event_params['pix_center'][0]
    bkg_events["ypix"] = pix[1, :] + event_params['pix_center'][1]

    mylog.info("Scattering energies with RMF %s." % os.path.split(rmf.filename)[-1])
    bkg_events = rmf.scatter_energies(bkg_events, prng=prng)

    return bkg_events