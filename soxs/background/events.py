import os

import numpy as np
from astropy.io import fits

from soxs.utils import get_rot_mat, mylog

key_map = {
    "telescope": "TELESCOP",
    "mission": "MISSION",
    "instrument": "INSTRUME",
    "channel_type": "CHANTYPE",
    "nchan": "PHA_BINS",
}


def add_background_from_file(events, event_params, bkg_file):
    from soxs.instrument import perform_dither

    f = fits.open(bkg_file)

    hdu = f["EVENTS"]

    dither_params = {}
    if "DITHXAMP" in hdu.header:
        dither_params["x_amp"] = hdu.header["DITHXAMP"]
        dither_params["y_amp"] = hdu.header["DITHYAMP"]
        dither_params["x_period"] = hdu.header["DITHXPER"]
        dither_params["y_period"] = hdu.header["DITHYPER"]
        dither_params["plate_scale"] = hdu.header["TCDLT3"] * 3600.0
        dither_params["dither_on"] = True
    else:
        dither_params["dither_on"] = False

    sexp = event_params["exposure_time"]
    bexp = hdu.header["EXPOSURE"]

    if event_params["exposure_time"] > hdu.header["EXPOSURE"]:
        raise RuntimeError(
            f"The background file does not have sufficient "
            f"exposure! Source exposure time {sexp}, background "
            f" exposure time {bexp}."
        )

    for k1, k2 in key_map.items():
        if event_params[k1] != hdu.header[k2]:
            raise RuntimeError(
                f"'{k1}' keyword does not match! "
                f"{event_params[k1]} vs. {hdu.header[k2]}"
            )
    rmf1 = os.path.split(event_params["rmf"])[-1]
    rmf2 = hdu.header["RESPFILE"]
    arf1 = os.path.split(event_params["arf"])[-1]
    arf2 = hdu.header["ANCRFILE"]
    if rmf1 != rmf2:
        raise RuntimeError(f"RMFs do not match! {rmf1} vs. {rmf2}")
    if arf1 != arf2:
        raise RuntimeError(f"ARFs do not match! {arf1} vs. {arf2}")

    idxs = hdu.data["TIME"] < sexp

    mylog.info("Adding %s background events from %s.", idxs.sum(), bkg_file)

    if event_params["roll_angle"] == hdu.header["ROLL_PNT"]:
        xpix = hdu.data["X"][idxs]
        ypix = hdu.data["Y"][idxs]
    else:
        rot_mat = get_rot_mat(event_params["roll_angle"])
        if dither_params["dither_on"]:
            t = hdu.data["TIME"][idxs]
            x_off, y_off = perform_dither(t, dither_params)
        else:
            x_off = 0.0
            y_off = 0.0
        det = np.array(
            [
                hdu.data["DETX"][idxs]
                + x_off
                - event_params["aimpt_coords"][0]
                - event_params["aimpt_shift"][0],
                hdu.data["DETY"][idxs]
                + y_off
                - event_params["aimpt_coords"][1]
                - event_params["aimpt_shift"][1],
            ]
        )
        xpix, ypix = np.dot(rot_mat.T, det)

        xpix += hdu.header["TCRPX2"]
        ypix += hdu.header["TCRPX3"]

    all_events = {}
    for key in ["detx", "dety", "time", "ccd_id", event_params["channel_type"]]:
        all_events[key] = np.concatenate([events[key], hdu.data[key.upper()][idxs]])
    all_events["xpix"] = np.concatenate([events["xpix"], xpix])
    all_events["ypix"] = np.concatenate([events["ypix"], ypix])
    all_events["energy"] = np.concatenate(
        [events["energy"], hdu.data["ENERGY"][idxs] * 1.0e-3]
    )

    f.close()

    return all_events
