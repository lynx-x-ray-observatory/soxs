import numpy as np
from astropy import wcs
from astropy.io import fits
from pathlib import Path, PurePath


def wcs_from_header(h):
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [h["TCRVL2"], h["TCRVL3"]]
    w.wcs.crpix = [h["TCRPX2"], h["TCRPX3"]]
    w.wcs.cdelt = [h["TCDLT2"], h["TCDLT3"]]
    w.wcs.ctype = [h["TCTYP2"], h["TCTYP3"]]
    w.wcs.cunit = [h["TCUNI2"], h["TCUNI3"]]
    return w


def make_event_file(events, parameters):
    from astropy.time import Time, TimeDelta

    t_begin = Time.now()
    dt = TimeDelta(parameters["exposure_time"], format="sec")
    t_end = t_begin + dt

    col_x = fits.Column(name="X", format="D", unit="pixel", array=events["xpix"])
    col_y = fits.Column(name="Y", format="D", unit="pixel", array=events["ypix"])
    col_e = fits.Column(
        name="ENERGY", format="E", unit="eV", array=events["energy"] * 1000.0
    )
    col_dx = fits.Column(name="DETX", format="D", unit="pixel", array=events["detx"])
    col_dy = fits.Column(name="DETY", format="D", unit="pixel", array=events["dety"])
    col_id = fits.Column(
        name="CCD_ID", format="J", unit="pixel", array=events["ccd_id"]
    )
    col_se = fits.Column(
        name="SOXS_ENERGY", format="E", unit="eV", array=events["soxs_energy"] * 1000.0
    )

    chantype = parameters["channel_type"].upper()
    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"
    col_ch = fits.Column(name=chantype, format="1J", unit=cunit, array=events[chantype])

    col_t = fits.Column(name="TIME", format="1D", unit="s", array=events["time"])

    cols = [col_e, col_x, col_y, col_ch, col_t, col_dx, col_dy, col_id, col_se]

    coldefs = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "EVENTS"

    tbhdu.header["MTYPE1"] = "sky"
    tbhdu.header["MFORM1"] = "x,y"
    tbhdu.header["MTYPE2"] = "EQPOS"
    tbhdu.header["MFORM2"] = "RA,DEC"
    tbhdu.header["TCTYP2"] = "RA---TAN"
    tbhdu.header["TCTYP3"] = "DEC--TAN"
    tbhdu.header["TCRVL2"] = parameters["sky_center"][0]
    tbhdu.header["TCRVL3"] = parameters["sky_center"][1]
    tbhdu.header["TCDLT2"] = -parameters["plate_scale"]
    tbhdu.header["TCDLT3"] = parameters["plate_scale"]
    tbhdu.header["TCRPX2"] = parameters["pix_center"][0]
    tbhdu.header["TCRPX3"] = parameters["pix_center"][1]
    tbhdu.header["TCUNI2"] = "deg"
    tbhdu.header["TCUNI3"] = "deg"
    tbhdu.header["TLMIN2"] = 0.5
    tbhdu.header["TLMIN3"] = 0.5
    tbhdu.header["TLMAX2"] = 2.0 * parameters["num_pixels"] + 0.5
    tbhdu.header["TLMAX3"] = 2.0 * parameters["num_pixels"] + 0.5
    tbhdu.header["TLMIN4"] = parameters["chan_lim"][0]
    tbhdu.header["TLMAX4"] = parameters["chan_lim"][1]
    tbhdu.header["TLMIN6"] = -0.5 * parameters["num_pixels"]
    tbhdu.header["TLMAX6"] = 0.5 * parameters["num_pixels"]
    tbhdu.header["TLMIN7"] = -0.5 * parameters["num_pixels"]
    tbhdu.header["TLMAX7"] = 0.5 * parameters["num_pixels"]
    tbhdu.header["EXPOSURE"] = parameters["exposure_time"]
    tbhdu.header["TSTART"] = 0.0
    tbhdu.header["TSTOP"] = parameters["exposure_time"]
    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["RADECSYS"] = "FK5"
    tbhdu.header["EQUINOX"] = 2000.0
    tbhdu.header["HDUCLASS"] = "OGIP"
    tbhdu.header["HDUCLAS1"] = "EVENTS"
    tbhdu.header["HDUCLAS2"] = "ACCEPTED"
    tbhdu.header["DATE"] = t_begin.tt.isot
    tbhdu.header["DATE-OBS"] = t_begin.tt.isot
    tbhdu.header["DATE-END"] = t_end.tt.isot
    tbhdu.header["RESPFILE"] = PurePath(parameters["rmf"]).parts[-1]
    tbhdu.header["PHA_BINS"] = parameters["nchan"]
    tbhdu.header["ANCRFILE"] = PurePath(parameters["arf"]).parts[-1]
    tbhdu.header["CHANTYPE"] = parameters["channel_type"]
    tbhdu.header["MISSION"] = parameters["mission"]
    tbhdu.header["TELESCOP"] = parameters["telescope"]
    tbhdu.header["INSTRUME"] = parameters["instrument"]
    tbhdu.header["RA_PNT"] = parameters["sky_center"][0]
    tbhdu.header["DEC_PNT"] = parameters["sky_center"][1]
    tbhdu.header["ROLL_PNT"] = parameters["roll_angle"]
    tbhdu.header["AIMPT_X"] = parameters["aimpt_coords"][0]
    tbhdu.header["AIMPT_Y"] = parameters["aimpt_coords"][1]
    tbhdu.header["AIMPT_DX"] = parameters["aimpt_shift"][0]
    tbhdu.header["AIMPT_DY"] = parameters["aimpt_shift"][1]
    if parameters["dither_params"]["dither_on"]:
        tbhdu.header["DITHXAMP"] = parameters["dither_params"]["x_amp"]
        tbhdu.header["DITHYAMP"] = parameters["dither_params"]["y_amp"]
        tbhdu.header["DITHXPER"] = parameters["dither_params"]["x_period"]
        tbhdu.header["DITHYPER"] = parameters["dither_params"]["y_period"]

    start = fits.Column(name="START", format="1D", unit="s", array=np.array([0.0]))
    stop = fits.Column(
        name="STOP",
        format="1D",
        unit="s",
        array=np.array([parameters["exposure_time"]]),
    )

    tbhdu_gti = fits.BinTableHDU.from_columns([start, stop])
    tbhdu_gti.name = "STDGTI"
    tbhdu_gti.header["TSTART"] = 0.0
    tbhdu_gti.header["TSTOP"] = parameters["exposure_time"]
    tbhdu_gti.header["HDUCLASS"] = "OGIP"
    tbhdu_gti.header["HDUCLAS1"] = "GTI"
    tbhdu_gti.header["HDUCLAS2"] = "STANDARD"
    tbhdu_gti.header["RADECSYS"] = "FK5"
    tbhdu_gti.header["EQUINOX"] = 2000.0
    tbhdu_gti.header["DATE"] = t_begin.tt.isot
    tbhdu_gti.header["DATE-OBS"] = t_begin.tt.isot
    tbhdu_gti.header["DATE-END"] = t_end.tt.isot

    hdulist = [fits.PrimaryHDU(), tbhdu, tbhdu_gti]

    return fits.HDUList(hdulist)


def _region_filter(hdu, region, format="ds9", exclude=False):
    from regions import PixCoord, PixelRegion, Region, Regions, SkyRegion

    if isinstance(region, str):
        if Path(region).exists():
            region = Regions.read(region, format=format)
        else:
            region = Regions.parse(region, format=format)
    elif not isinstance(region, (Region, Regions)):
        raise RuntimeError("'region' argument is not valid!")
    pixcoords = PixCoord(hdu.data["X"], hdu.data["Y"])
    if isinstance(region, Region):
        region = [region]
    evt_mask = np.zeros(hdu.data["ENERGY"].size, dtype="bool")
    for r in region:
        include_this = bool(r.meta.get("include", True))
        if isinstance(r, PixelRegion):
            this_mask = r.contains(pixcoords)
        elif isinstance(r, SkyRegion):
            w = wcs_from_header(hdu.header)
            skycoords = pixcoords.to_sky(w, origin=1)
            this_mask = r.contains(skycoords, w)
        else:
            raise NotImplementedError
        if include_this:
            evt_mask |= this_mask
        else:
            evt_mask &= this_mask
    if exclude:
        evt_mask = ~evt_mask
    return evt_mask


def _combine_events(eventfiles, wcs_out, shape_out, outfile, overwrite=False):
    x = []
    y = []
    e = []
    ch = []
    t = []
    se = []
    chantype = ""
    header_dict = {}
    for i, eventfile in enumerate(eventfiles):
        with fits.open(eventfile, memmap=True) as f:
            hdu = f["EVENTS"]
            wcs_in = wcs_from_header(hdu.header)
            if i == 0:
                chantype = hdu.header["CHANTYPE"]
                for key in [
                    "TELESCOP",
                    "INSTRUME",
                    "ANCRFILE",
                    "RESPFILE",
                    "EXPOSURE",
                    "CHANTYPE",
                    "DATE",
                    "DATE-OBS",
                    "DATE-END",
                    "MISSION",
                    "TSTART",
                    "TSTOP",
                    "TLMIN4",
                    "TLMAX4",
                    "PHA_BINS",
                ]:
                    header_dict[key] = hdu.header[key]
                idxs = ...
            else:
                for key in [
                    "TELESCOP",
                    "INSTRUME",
                    "ANCRFILE",
                    "RESPFILE",
                    "MISSION",
                    "TLMIN4",
                    "TLMAX4",
                    "PHA_BINS",
                    "CHANTYPE",
                ]:
                    v1 = hdu.header[key]
                    v2 = header_dict[key]
                    if v1 != v2:
                        raise ValueError(
                            f"'{key}' from {eventfile} "
                            f"does not match from {eventfiles[0]}! "
                            f"{v1} != {v2}!"
                        )
                if hdu.header["EXPOSURE"] < header_dict["EXPOSURE"]:
                    raise ValueError(
                        f"Exposure time for {eventfile} is less than "
                        f"that for {eventfiles[0]}! {hdu.header['EXPOSURE']} "
                        f"< {header_dict['EXPOSURE']}!"
                    )
                idxs = hdu.data["TIME"] < header_dict["EXPOSURE"]

            ra, dec = wcs_in.wcs_pix2world(hdu.data["X"][idxs], hdu.data["Y"][idxs], 1)
            xx, yy = wcs_out.wcs_world2pix(ra, dec, 1)
            x.append(xx)
            y.append(yy)
            e.append(hdu.data["ENERGY"][idxs])
            if "SOXS_ENERGY" in hdu.data.dtype.names:
                se.append(hdu.data["SOXS_ENERGY"][idxs])
            ch.append(hdu.data[hdu.header["CHANTYPE"]][idxs])
            t.append(hdu.data["TIME"][idxs])

    x = np.concatenate(x)
    y = np.concatenate(y)
    e = np.concatenate(e)
    se = np.concatenate(se)
    ch = np.concatenate(ch)
    t = np.concatenate(t)

    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"

    col_x = fits.Column(name="X", format="D", unit="pixel", array=x)
    col_y = fits.Column(name="Y", format="D", unit="pixel", array=y)
    col_e = fits.Column(name="ENERGY", format="E", unit="eV", array=e)
    col_ch = fits.Column(name=chantype, format="1J", unit=cunit, array=ch)
    col_t = fits.Column(name="TIME", format="1D", unit="s", array=t)
    coldefs = [col_e, col_x, col_y, col_ch, col_t]
    if len(se) == len(e):
        col_se = fits.Column(name="SOXS_ENERGY", format="E", unit="eV", array=se)
        coldefs.append(col_se)
    coldefs = fits.ColDefs(coldefs)
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "EVENTS"

    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["RADECSYS"] = "FK5"
    tbhdu.header["EQUINOX"] = 2000.0
    tbhdu.header["HDUCLASS"] = "OGIP"
    tbhdu.header["HDUCLAS1"] = "EVENTS"
    tbhdu.header["HDUCLAS2"] = "ACCEPTED"
    tbhdu.header["MTYPE1"] = "sky"
    tbhdu.header["MFORM1"] = "x,y"
    tbhdu.header["MTYPE2"] = "EQPOS"
    tbhdu.header["MFORM2"] = "RA,DEC"
    tbhdu.header["TCTYP2"] = "RA---TAN"
    tbhdu.header["TCTYP3"] = "DEC--TAN"
    tbhdu.header["TCRVL2"] = wcs_out.wcs.crval[0]
    tbhdu.header["TCRVL3"] = wcs_out.wcs.crval[1]
    tbhdu.header["TCDLT2"] = wcs_out.wcs.cdelt[0]
    tbhdu.header["TCDLT3"] = wcs_out.wcs.cdelt[1]
    tbhdu.header["TCRPX2"] = wcs_out.wcs.crpix[0]
    tbhdu.header["TCRPX3"] = wcs_out.wcs.crpix[1]
    tbhdu.header["TCUNI2"] = "deg"
    tbhdu.header["TCUNI3"] = "deg"
    tbhdu.header["TLMIN2"] = 0.5
    tbhdu.header["TLMIN3"] = 0.5
    tbhdu.header["TLMAX2"] = shape_out[0] + 0.5
    tbhdu.header["TLMAX3"] = shape_out[1] + 0.5
    tbhdu.header["CHANTYPE"] = chantype

    with fits.open(eventfiles[0], memmap=True) as f:
        for key in [
            "TELESCOP",
            "INSTRUME",
            "ANCRFILE",
            "RESPFILE",
            "EXPOSURE",
            "DATE",
            "DATE-OBS",
            "DATE-END",
            "MISSION",
            "TSTART",
            "TSTOP",
            "TLMIN4",
            "TLMAX4",
            "PHA_BINS",
        ]:
            tbhdu.header[key] = f["EVENTS"].header[key]

    start = fits.Column(name="START", format="1D", unit="s", array=np.array([0.0]))
    stop = fits.Column(
        name="STOP", format="1D", unit="s", array=np.array([tbhdu.header["EXPOSURE"]])
    )

    tbhdu_gti = fits.BinTableHDU.from_columns([start, stop])
    tbhdu_gti.name = "STDGTI"
    tbhdu_gti.header["TSTART"] = 0.0
    tbhdu_gti.header["TSTOP"] = tbhdu.header["EXPOSURE"]
    tbhdu_gti.header["HDUCLASS"] = "OGIP"
    tbhdu_gti.header["HDUCLAS1"] = "GTI"
    tbhdu_gti.header["HDUCLAS2"] = "STANDARD"
    tbhdu_gti.header["RADECSYS"] = "FK5"
    tbhdu_gti.header["EQUINOX"] = 2000.0
    tbhdu_gti.header["DATE"] = tbhdu.header["DATE"]
    tbhdu_gti.header["DATE-OBS"] = tbhdu.header["DATE-OBS"]
    tbhdu_gti.header["DATE-END"] = tbhdu.header["DATE-END"]

    hdulist = [fits.PrimaryHDU(), tbhdu, tbhdu_gti]

    fits.HDUList(hdulist).writeto(outfile, overwrite=overwrite)
