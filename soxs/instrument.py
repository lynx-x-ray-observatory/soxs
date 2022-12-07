import os
import warnings
from collections import defaultdict

import numpy as np
from astropy import wcs
from regions import PixCoord

from soxs.events import write_event_file
from soxs.instrument_registry import instrument_registry
from soxs.psf import psf_model_registry
from soxs.response import AuxiliaryResponseFile, RedistributionMatrixFile
from soxs.simput import SimputPhotonList, read_simput_catalog
from soxs.utils import (
    create_region,
    ensure_numpy_array,
    get_data_file,
    get_rot_mat,
    mylog,
    parse_prng,
    parse_value,
)


def perform_dither(t, dither_dict):
    if dither_dict["dither_on"]:
        a = 2.0 * np.pi / dither_dict["x_period"]
        b = 2.0 * np.pi / dither_dict["y_period"]
        A = dither_dict["x_amp"] / dither_dict["plate_scale"]
        B = dither_dict["y_amp"] / dither_dict["plate_scale"]
        x_offset = A * np.sin(a * t)
        y_offset = B * np.sin(b * t)
    else:
        x_offset = np.zeros(t.size)
        y_offset = np.zeros(t.size)
    return x_offset, y_offset


def make_source_list(source):
    if source is None:
        source_list = []
        parameters = {}
    elif isinstance(source, dict):
        parameters = {}
        for key in ["flux", "emin", "emax", "src_names"]:
            parameters[key] = source[key]
        source_list = []
        for i in range(len(parameters["flux"])):
            phlist = SimputPhotonList(
                source["ra"][i],
                source["dec"][i],
                source["energy"][i],
                parameters["flux"][i],
                parameters["src_names"][i],
            )
            source_list.append(phlist)
    elif isinstance(source, str):
        # Assume this is a SIMPUT catalog
        source_list, parameters = read_simput_catalog(source)
    return source_list, parameters


def generate_events(
    source,
    exp_time,
    instrument,
    sky_center,
    no_dither=False,
    dither_params=None,
    roll_angle=0.0,
    subpixel_res=False,
    aimpt_shift=None,
    prng=None,
):
    """
    Take unconvolved events and convolve them with instrumental responses. This
    function does the following:

    1. Determines which events are observed using the ARF
    2. Pixelizes the events, applying PSF effects and dithering
    3. Determines energy channels using the RMF

    This function is not meant to be called by the end-user but is used by
    the :func:`~soxs.instrument.instrument_simulator` function.

    Parameters
    ----------
    input_events : string, dict, or None
        The unconvolved events to be used as input. Can be one of the
        following:
        1. The name of a SIMPUT catalog file.
        2. A Python dictionary containing the following items:
        "ra": A NumPy array of right ascension values in degrees.
        "dec": A NumPy array of declination values in degrees.
        "energy": A NumPy array of energy values in keV.
        "flux": The flux of the entire source, in units of erg/cm**2/s.
    out_file : string
        The name of the event file to be written.
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time to use, in seconds.
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry.
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    no_dither : boolean, optional
        If True, turn off dithering entirely. Default: False
    dither_params : array-like of floats, optional
        The parameters to use to control the size and period of the dither
        pattern. The first two numbers are the dither amplitude in x and y
        detector coordinates in arcseconds, and the second two numbers are
        the dither period in x and y detector coordinates in seconds.
        Default: [8.0, 8.0, 1000.0, 707.0].
    roll_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The roll angle of the observation in degrees. Default: 0.0
    subpixel_res : boolean, optional
        If True, event positions are not randomized within the pixels
        within which they are detected. Default: False
    aimpt_shift : array-like, optional
        A two-float array-like object which shifts the aimpoint on the
        detector from the nominal position. Units are in arcseconds.
        Default: None, which results in no shift from the nominal aimpoint.
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    """
    exp_time = parse_value(exp_time, "s")
    roll_angle = parse_value(roll_angle, "deg")
    prng = parse_prng(prng)
    source_list, parameters = make_source_list(source)

    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError(f"Instrument {instrument} is not in the instrument registry!")
    if not instrument_spec["imaging"]:
        raise RuntimeError(
            f"Instrument '{instrument_spec['name']}' is not "
            f"designed for imaging observations!"
        )

    arf_file = get_data_file(instrument_spec["arf"])
    rmf_file = get_data_file(instrument_spec["rmf"])
    arf = AuxiliaryResponseFile(arf_file)
    rmf = RedistributionMatrixFile(rmf_file)

    nx = instrument_spec["num_pixels"]
    plate_scale = instrument_spec["fov"] / nx / 60.0  # arcmin to deg
    plate_scale_arcsec = plate_scale * 3600.0

    if aimpt_shift is None:
        aimpt_shift = np.zeros(2)
    aimpt_shift = ensure_numpy_array(aimpt_shift).astype("float64")
    aimpt_shift /= plate_scale_arcsec

    if not instrument_spec["dither"]:
        dither_on = False
    else:
        dither_on = not no_dither
    if dither_params is None:
        dither_params = [8.0, 8.0, 1000.0, 707.0]
    dither_dict = {
        "x_amp": dither_params[0],
        "y_amp": dither_params[1],
        "x_period": dither_params[2],
        "y_period": dither_params[3],
        "dither_on": dither_on,
        "plate_scale": plate_scale_arcsec,
    }

    event_params = {
        "exposure_time": exp_time,
        "arf": arf.filename,
        "sky_center": sky_center,
        "pix_center": np.array([0.5 * (2 * nx + 1)] * 2),
        "num_pixels": nx,
        "plate_scale": plate_scale,
        "rmf": rmf.filename,
        "channel_type": rmf.chan_type.upper(),
        "telescope": rmf.header["TELESCOP"],
        "instrument": instrument_spec["name"],
        "mission": rmf.header.get("MISSION", ""),
        "nchan": rmf.n_ch,
        "roll_angle": roll_angle,
        "fov": instrument_spec["fov"],
        "chan_lim": [rmf.cmin, rmf.cmax],
        "chips": instrument_spec["chips"],
        "dither_params": dither_dict,
        "aimpt_coords": instrument_spec["aimpt_coords"],
        "aimpt_shift": aimpt_shift,
    }

    # Set up WCS

    w = wcs.WCS(naxis=2)
    w.wcs.crval = event_params["sky_center"]
    w.wcs.crpix = event_params["pix_center"]
    w.wcs.cdelt = [-plate_scale, plate_scale]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ["deg"] * 2

    # Determine rotation matrix
    rot_mat = get_rot_mat(roll_angle)

    # Set up PSF
    psf_type = instrument_spec["psf"][0]
    psf_class = psf_model_registry[psf_type]
    psf = psf_class(instrument_spec, prng=prng)

    all_events = defaultdict(list)

    for i, src in enumerate(source_list):

        mylog.info("Detecting events from source %s.", parameters["src_names"][i])

        # Step 1: Use ARF to determine which photons are observed

        mylog.info(
            "Applying energy-dependent effective area from %s.",
            os.path.split(arf.filename)[-1],
        )
        refband = [parameters["emin"][i], parameters["emax"][i]]
        if src.src_type == "phlist":
            events = arf.detect_events_phlist(
                src.events.copy(), exp_time, parameters["flux"][i], refband, prng=prng
            )
        elif src.src_type.endswith("spectrum"):
            events = arf.detect_events_spec(src, exp_time, refband, prng=prng)

        n_evt = events["energy"].size

        if n_evt == 0:
            mylog.warning("No events were observed for this source!!!")
        else:

            # Step 2: Assign pixel coordinates to events. Apply dithering and
            # PSF. Clip events that don't fall within the detection region.

            mylog.info("Pixeling events.")

            # Convert RA, Dec to pixel coordinates
            xpix, ypix = w.wcs_world2pix(events["ra"], events["dec"], 1)

            xpix -= event_params["pix_center"][0]
            ypix -= event_params["pix_center"][1]

            events.pop("ra")
            events.pop("dec")

            n_evt = xpix.size

            # Rotate physical coordinates to detector coordinates

            det = np.dot(rot_mat, np.array([xpix, ypix]))
            detx = det[0, :] + event_params["aimpt_coords"][0] + aimpt_shift[0]
            dety = det[1, :] + event_params["aimpt_coords"][1] + aimpt_shift[1]

            # Add times to events
            events["time"] = prng.uniform(
                size=n_evt, low=0.0, high=event_params["exposure_time"]
            )

            # Apply dithering

            x_offset, y_offset = perform_dither(events["time"], dither_dict)

            detx -= x_offset
            dety -= y_offset

            # PSF scattering of detector coordinates

            mylog.info("Scattering events with a %s-based PSF.", psf)
            detx, dety = psf.scatter(detx, dety, events["energy"])

            # Convert detector coordinates to chip coordinates.
            # Throw out events that don't fall on any chip.

            cx = np.trunc(detx) + 0.5 * np.sign(detx)
            cy = np.trunc(dety) + 0.5 * np.sign(dety)

            events["chip_id"] = -np.ones(n_evt, dtype="int")
            for i, chip in enumerate(event_params["chips"]):
                rtype = chip[0]
                args = chip[1:]
                r, _ = create_region(rtype, args, 0.0, 0.0)
                inside = r.contains(PixCoord(cx, cy))
                events["chip_id"][inside] = i
            keep = events["chip_id"] > -1

            mylog.debug(
                "%d events were rejected because they do not fall on any CCD.",
                n_evt - keep.sum(),
            )
            n_evt = keep.sum()

            if n_evt == 0:
                mylog.warning(
                    "No events are within the field " "of view for this source!!!"
                )
            else:

                mylog.info("%d events were detected from the source.", n_evt)

                # Keep only those events which fall on a chip

                for key in events:
                    events[key] = events[key][keep]

                # Convert chip coordinates back to detector coordinates,
                # unless the user has specified that they want subpixel
                # resolution

                if subpixel_res:
                    events["detx"] = detx[keep]
                    events["dety"] = dety[keep]
                else:
                    events["detx"] = cx[keep] + prng.uniform(
                        low=-0.5, high=0.5, size=n_evt
                    )
                    events["dety"] = cy[keep] + prng.uniform(
                        low=-0.5, high=0.5, size=n_evt
                    )

                # Convert detector coordinates back to pixel coordinates by
                # adding the dither offsets back in and applying the rotation
                # matrix again

                det = np.array(
                    [
                        events["detx"]
                        + x_offset[keep]
                        - event_params["aimpt_coords"][0]
                        - aimpt_shift[0],
                        events["dety"]
                        + y_offset[keep]
                        - event_params["aimpt_coords"][1]
                        - aimpt_shift[1],
                    ]
                )
                pix = np.dot(rot_mat.T, det)

                events["xpix"] = pix[0, :] + event_params["pix_center"][0]
                events["ypix"] = pix[1, :] + event_params["pix_center"][1]

        if n_evt > 0:
            for key in events:
                all_events[key] = np.concatenate([all_events[key], events[key]])

    if len(all_events["energy"]) == 0:
        mylog.warning(
            "No events from any of the sources in " "the catalog were detected!"
        )
        for key in [
            "xpix",
            "ypix",
            "detx",
            "dety",
            "time",
            "chip_id",
            event_params["channel_type"],
        ]:
            all_events[key] = np.array([])
    else:
        # Step 4: Scatter energies with RMF
        mylog.info("Scattering energies with RMF %s.", os.path.split(rmf.filename)[-1])
        all_events = rmf.scatter_energies(all_events, prng=prng)

    return all_events, event_params


def make_background(
    exp_time,
    instrument,
    sky_center,
    foreground=True,
    ptsrc_bkgnd=True,
    instr_bkgnd=True,
    no_dither=False,
    dither_params=None,
    roll_angle=0.0,
    subpixel_res=False,
    input_pt_sources=None,
    aimpt_shift=None,
    prng=None,
    **kwargs,
):
    """
    Make background events.

    Parameters
    ----------
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time to use, in seconds.
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry.
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    foreground : boolean, optional
        Whether to include the Galactic foreground. Default: True
    instr_bkgnd : boolean, optional
        Whether to include the instrumental background. Default: True
    no_dither : boolean, optional
        If True, turn off dithering entirely. Default: False
    dither_params : array-like of floats, optional
        The parameters to use to control the size and period of the dither
        pattern. The first two numbers are the dither amplitude in x and y
        detector coordinates in arcseconds, and the second two numbers are
        the dither period in x and y detector coordinates in seconds.
        Default: [8.0, 8.0, 1000.0, 707.0].
    ptsrc_bkgnd : boolean, optional
        Whether to include the point-source background. Default: True
    roll_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The roll angle of the observation in degrees. Default: 0.0
    subpixel_res: boolean, optional
        If True, event positions are not randomized within the pixels
        within which they are detected. Default: False
    input_pt_sources : string, optional
        If set to a filename, input the point source positions, fluxes,
        and spectral indices from an ASCII table instead of generating
        them. Default: None
    aimpt_shift : array-like, optional
        A two-float array-like object which shifts the aimpoint on the
        detector from the nominal position. Units are in arcseconds.
        Default: None, which results in no shift from the nominal aimpoint.
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically this will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    """
    if "nH" in kwargs or "absorb_model" in kwargs:
        warnings.warn(
            "The 'nH' and 'absorb_model' keyword arguments"
            "have been omitted. Please set the 'bkgnd_nH' "
            "and 'bkgnd_absorb_model' values in the SOXS"
            "configuration file if you want to change these "
            "values. ",
            DeprecationWarning,
        )
    from soxs.background import make_diffuse_background, make_ptsrc_background

    prng = parse_prng(prng)
    exp_time = parse_value(exp_time, "s")
    roll_angle = parse_value(roll_angle, "deg")
    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError(
            f"Instrument {instrument} is not in the " f"instrument registry!"
        )
    if not instrument_spec["imaging"]:
        raise RuntimeError(
            f"Instrument '{instrument_spec['name']}' is not "
            f"designed for imaging observations!"
        )
    fov = instrument_spec["fov"]

    input_events = defaultdict(list)

    arf_file = get_data_file(instrument_spec["arf"])
    arf = AuxiliaryResponseFile(arf_file)
    rmf_file = get_data_file(instrument_spec["rmf"])
    rmf = RedistributionMatrixFile(rmf_file)

    if ptsrc_bkgnd:
        mylog.info("Adding in point-source background.")
        ptsrc_events = make_ptsrc_background(
            exp_time,
            fov,
            sky_center,
            area=1.2 * arf.max_area,
            input_sources=input_pt_sources,
            prng=prng,
        )
        for key in ["ra", "dec", "energy"]:
            input_events[key].append(ptsrc_events[key])
        input_events["flux"].append(ptsrc_events["flux"])
        input_events["emin"].append(ptsrc_events["energy"].min())
        input_events["emax"].append(ptsrc_events["energy"].max())
        input_events["src_names"].append("ptsrc_bkgnd")
        events, event_params = generate_events(
            input_events,
            exp_time,
            instrument,
            sky_center,
            no_dither=no_dither,
            dither_params=dither_params,
            roll_angle=roll_angle,
            subpixel_res=subpixel_res,
            aimpt_shift=aimpt_shift,
            prng=prng,
        )
        mylog.info(
            "Generated %d photons from the point-source background.",
            events["energy"].size,
        )
    else:
        nx = instrument_spec["num_pixels"]
        plate_scale = instrument_spec["fov"] / nx / 60.0
        plate_scale_arcsec = plate_scale * 3600.0
        if aimpt_shift is None:
            aimpt_shift = np.zeros(2)
        aimpt_shift = ensure_numpy_array(aimpt_shift).astype("float64")
        aimpt_shift /= plate_scale_arcsec
        events = defaultdict(list)
        if not instrument_spec["dither"]:
            dither_on = False
        else:
            dither_on = not no_dither
        if dither_params is None:
            dither_params = [8.0, 8.0, 1000.0, 707.0]
        dither_dict = {
            "x_amp": dither_params[0],
            "y_amp": dither_params[1],
            "x_period": dither_params[2],
            "y_period": dither_params[3],
            "dither_on": dither_on,
            "plate_scale": instrument_spec["fov"] / nx * 60.0,
        }
        event_params = {
            "exposure_time": exp_time,
            "fov": instrument_spec["fov"],
            "num_pixels": nx,
            "pix_center": np.array([0.5 * (2 * nx + 1)] * 2),
            "channel_type": rmf.header["CHANTYPE"].upper(),
            "sky_center": sky_center,
            "dither_params": dither_dict,
            "plate_scale": plate_scale,
            "chan_lim": [rmf.cmin, rmf.cmax],
            "rmf": rmf_file,
            "arf": arf_file,
            "telescope": rmf.header["TELESCOP"],
            "instrument": instrument_spec["name"],
            "mission": rmf.header.get("MISSION", ""),
            "nchan": rmf.n_ch,
            "roll_angle": roll_angle,
            "aimpt_coords": instrument_spec["aimpt_coords"],
            "aimpt_shift": aimpt_shift,
        }

    if "chips" not in event_params:
        event_params["chips"] = instrument_spec["chips"]

    instr_bkgnd &= instrument_spec["bkgnd"] is not None

    if foreground or instr_bkgnd:
        bkg_events = make_diffuse_background(
            foreground, instr_bkgnd, instrument_spec, event_params, arf, rmf, prng=prng
        )
        for key in bkg_events:
            events[key] = np.concatenate([events[key], bkg_events[key]])

    return events, event_params


def make_background_file(
    out_file,
    exp_time,
    instrument,
    sky_center,
    overwrite=False,
    foreground=True,
    instr_bkgnd=True,
    ptsrc_bkgnd=True,
    no_dither=False,
    dither_params=None,
    subpixel_res=False,
    input_pt_sources=None,
    prng=None,
    **kwargs,
):
    """
    Make an event file consisting entirely of background events. This will be
    useful for creating backgrounds that can be added to simulations of sources.

    Parameters
    ----------
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time to use, in seconds.
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry.
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    overwrite : boolean, optional
        Whether to overwrite an existing file with the same name.
        Default: False
    foreground : boolean, optional
        Whether or not to include the Galactic foreground. Default: True
    instr_bkgnd : boolean, optional
        Whether or not to include the instrumental background. Default: True
    ptsrc_bkgnd : boolean, optional
        Whether or not to include the point-source background. Default: True
    no_dither : boolean, optional
        If True, turn off dithering entirely. Default: False
    dither_params : array-like of floats, optional
        The parameters to use to control the size and period of the dither
        pattern. The first two numbers are the dither amplitude in x and y
        detector coordinates in arcseconds, and the second two numbers are
        the dither period in x and y detector coordinates in seconds.
        Default: [8.0, 8.0, 1000.0, 707.0].
    subpixel_res: boolean, optional
        If True, event positions are not randomized within the pixels
        within which they are detected. Default: False
    input_pt_sources : string, optional
        If set to a filename, input the point source positions, fluxes,
        and spectral indices from an ASCII table instead of generating
        them. Default: None
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.
    """
    if "nH" in kwargs or "absorb_model" in kwargs:
        warnings.warn(
            "The 'nH' and 'absorb_model' keyword arguments"
            "have been omitted. Please set the 'bkgnd_nH' "
            "and 'bkgnd_absorb_model' values in the SOXS"
            "configuration file if you want to change these "
            "values. ",
            DeprecationWarning,
        )
    if "input_sources" in kwargs:
        warnings.warn(
            "The 'input_sources' keyword argument has been changed "
            "to 'input_pt_sources' and is deprecated.",
            DeprecationWarning,
        )
        input_pt_sources = kwargs.pop("input_sources")
    prng = parse_prng(prng)
    events, event_params = make_background(
        exp_time,
        instrument,
        sky_center,
        ptsrc_bkgnd=ptsrc_bkgnd,
        foreground=foreground,
        instr_bkgnd=instr_bkgnd,
        no_dither=no_dither,
        dither_params=dither_params,
        subpixel_res=subpixel_res,
        input_pt_sources=input_pt_sources,
        prng=prng,
    )
    write_event_file(events, event_params, out_file, overwrite=overwrite)


def instrument_simulator(
    input_events,
    out_file,
    exp_time,
    instrument,
    sky_center,
    overwrite=False,
    instr_bkgnd=True,
    foreground=True,
    ptsrc_bkgnd=True,
    bkgnd_file=None,
    no_dither=False,
    dither_params=None,
    roll_angle=0.0,
    subpixel_res=False,
    aimpt_shift=None,
    input_pt_sources=None,
    prng=None,
):
    """
    Take unconvolved events and create an event file from them. This
    function calls generate_events to do the following:

    1. Determines which events are observed using the ARF
    2. Pixelizes the events, applying PSF effects and dithering
    3. Determines energy channels using the RMF

    and then calls make_background to add instrumental and astrophysical
    backgrounds, unless a background file is provided, in which case
    the background events are read from this file. The events are
    then written out to a file.

    Parameters
    ----------
    input_events : string, dict, or None
        The unconvolved events to be used as input. Can be one of the
        following:
        1. The name of a SIMPUT catalog file.
        2. A Python dictionary containing the following items:
        "ra": A NumPy array of right ascension values in degrees.
        "dec": A NumPy array of declination values in degrees.
        "energy": A NumPy array of energy values in keV.
        "flux": The flux of the entire source, in units of erg/cm**2/s.
    out_file : string
        The name of the event file to be written.
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time to use, in seconds.
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry.
    sky_center : array, tuple, or list
        The center RA, Dec coordinates of the observation, in degrees.
    overwrite : boolean, optional
        Whether to overwrite an existing file with the same name.
        Default: False
    instr_bkgnd : boolean, optional
        Whether to include the instrumental/particle background.
        Default: True
    foreground : boolean, optional
        Whether to include the local foreground.
        Default: True
    ptsrc_bkgnd : boolean, optional
        Whether to include the point-source background.
        Default: True
    bkgnd_file : string, optional
        If set, backgrounds will be loaded from this file and not generated
        on the fly. Default: None
    no_dither : boolean, optional
        If True, turn off dithering entirely. Default: False
    dither_params : array-like of floats, optional
        The parameters to use to control the size and period of the dither
        pattern. The first two numbers are the dither amplitude in x and y
        detector coordinates in arcseconds, and the second two numbers are
        the dither period in x and y detector coordinates in seconds.
        Default: [8.0, 8.0, 1000.0, 707.0].
    roll_angle : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The roll angle of the observation in degrees. Default: 0.0
    subpixel_res: boolean, optional
        If True, event positions are not randomized within the pixels
        within which they are detected. Default: False
    aimpt_shift : array-like, optional
        A two-float array-like object which shifts the aimpoint on the
        detector from the nominal position. Units are in arcseconds.
        Default: None, which results in no shift from the nominal aimpoint.
    input_pt_sources : string, optional
        If set to a filename, input the point source positions, fluxes,
        and spectral indices from an ASCII table instead of generating
        them. Default: None
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.

    Examples
    --------
    >>> instrument_simulator("sloshing_simput.fits", "sloshing_evt.fits",
    ...                      300000.0, "lynx_hdxi", [30., 45.], overwrite=True)
    """
    from soxs.background import add_background_from_file

    if not out_file.endswith(".fits"):
        out_file += ".fits"
    mylog.info("Making observation of source in %s.", out_file)
    # Make the source first
    events, event_params = generate_events(
        input_events,
        exp_time,
        instrument,
        sky_center,
        no_dither=no_dither,
        dither_params=dither_params,
        roll_angle=roll_angle,
        subpixel_res=subpixel_res,
        aimpt_shift=aimpt_shift,
        prng=prng,
    )
    # If the user wants backgrounds, either make the background or add an already existing
    # background event file. It may be necessary to reproject events to a new coordinate system.
    if bkgnd_file is None:
        if not instr_bkgnd and not ptsrc_bkgnd and not foreground:
            mylog.info("No backgrounds will be added to this observation.")
        else:
            mylog.info("Adding background events.")
            bkg_events, _ = make_background(
                exp_time,
                instrument,
                sky_center,
                foreground=foreground,
                instr_bkgnd=instr_bkgnd,
                no_dither=no_dither,
                dither_params=dither_params,
                ptsrc_bkgnd=ptsrc_bkgnd,
                prng=prng,
                subpixel_res=subpixel_res,
                roll_angle=roll_angle,
                aimpt_shift=aimpt_shift,
                input_pt_sources=input_pt_sources,
            )
            for key in events:
                events[key] = np.concatenate([events[key], bkg_events[key]])
    else:
        mylog.info("Adding background events from the file %s.", bkgnd_file)
        if not os.path.exists(bkgnd_file):
            raise IOError(f"Cannot find the background event file {bkgnd_file}!")
        events = add_background_from_file(events, event_params, bkgnd_file)
    if len(events["energy"]) == 0:
        mylog.warning(
            "No events were detected from source or background!! We "
            "will not write an event file."
        )
    else:
        write_event_file(events, event_params, out_file, overwrite=overwrite)
    mylog.info("Observation complete.")


def simulate_spectrum(
    spec,
    instrument,
    exp_time,
    out_file,
    instr_bkgnd=False,
    foreground=False,
    ptsrc_bkgnd=False,
    bkgnd_area=None,
    overwrite=False,
    prng=None,
    **kwargs,
):
    """
    Generate a PI or PHA spectrum from a :class:`~soxs.spectra.Spectrum`
    by convolving it with responses. To be used if one wants to
    create a spectrum without worrying about spatial response. Similar
    to XSPEC's "fakeit".

    Parameters
    ----------
    spec : :class:`~soxs.spectra.Spectrum`
        The spectrum to be convolved. If None is supplied, only backgrounds
        will be simulated (if they are turned on).
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry.
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time in seconds.
    out_file : string
        The file to write the spectrum to.
    instr_bkgnd : boolean, optional
        Whether to include the instrumental/particle background.
        Default: False
    foreground : boolean, optional
        Whether to include the local foreground.
        Default: False
    ptsrc_bkgnd : boolean, optional
        Whether to include the unresolved point-source background.
        Default: False
    bkgnd_area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The area on the sky for the background components, in square arcminutes.
        Default: None, necessary to specify if any of the background components
        are turned on.
    overwrite : boolean, optional
        Whether to overwrite an existing file. Default: False
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only
        be specified if you have a reason to generate the same
        set of random numbers, such as for a test. Default is None,
        which sets the seed based on the system time.

    Examples
    --------
    >>> spec = soxs.Spectrum.from_file("my_spectrum.txt")
    >>> soxs.simulate_spectrum(spec, "lynx_lxm", 100000.0,
    ...                        "my_spec.pi", overwrite=True)
    """
    from soxs.background.diffuse import (
        generate_channel_spectrum,
        make_frgnd_spectrum,
        read_instr_spectrum,
    )
    from soxs.background.spectra import BackgroundSpectrum
    from soxs.events import _write_spectrum
    from soxs.response import AuxiliaryResponseFile, RedistributionMatrixFile
    from soxs.spectra import ConvolvedSpectrum
    from soxs.utils import soxs_cfg

    if "nH" in kwargs or "absorb_model" in kwargs:
        warnings.warn(
            "The 'nH' and 'absorb_model' keyword arguments"
            "have been omitted. Please set the 'bkgnd_nH' "
            "and 'bkgnd_absorb_model' values in the SOXS"
            "configuration file if you want to change these "
            "values. ",
            DeprecationWarning,
        )
    prng = parse_prng(prng)
    exp_time = parse_value(exp_time, "s")
    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError(f"Instrument {instrument} is not in the instrument registry!")
    if foreground or instr_bkgnd or ptsrc_bkgnd:
        if instrument_spec["grating"]:
            raise NotImplementedError(
                "Backgrounds cannot be included in simulations "
                "of gratings spectra at this time!"
            )
        if bkgnd_area is None:
            raise RuntimeError(
                "The 'bkgnd_area' argument must be set if one wants "
                "to simulate backgrounds! Specify a value in square "
                "arcminutes."
            )
        bkgnd_area = np.sqrt(parse_value(bkgnd_area, "arcmin**2"))
    elif spec is None:
        raise RuntimeError("You have specified no source spectrum and no backgrounds!")
    arf_file = get_data_file(instrument_spec["arf"])
    rmf_file = get_data_file(instrument_spec["rmf"])
    arf = AuxiliaryResponseFile(arf_file)
    rmf = RedistributionMatrixFile(rmf_file)

    event_params = {
        "RESPFILE": os.path.split(rmf.filename)[-1],
        "ANCRFILE": os.path.split(arf.filename)[-1],
        "TELESCOP": rmf.header["TELESCOP"],
        "INSTRUME": rmf.header["INSTRUME"],
        "MISSION": rmf.header.get("MISSION", ""),
    }

    out_spec = np.zeros(rmf.n_ch)

    if spec is not None:
        cspec = ConvolvedSpectrum.convolve(spec, arf, use_arf_energies=True)
        out_spec += rmf.convolve_spectrum(cspec, exp_time, prng=prng)

    fov = None if bkgnd_area is None else np.sqrt(bkgnd_area)

    if foreground:
        mylog.info("Adding in astrophysical foreground.")
        frgnd_spec = rmf.convolve_spectrum(
            make_frgnd_spectrum(arf, rmf), exp_time, noisy=False, rate=True
        )
        out_spec += generate_channel_spectrum(
            frgnd_spec, exp_time, bkgnd_area, prng=prng
        )
    if instr_bkgnd and instrument_spec["bkgnd"] is not None:
        mylog.info("Adding in instrumental background.")
        bkgnd_spec = instrument_spec["bkgnd"]
        # Temporary hack for ACIS-S
        if "aciss" in instrument_spec["name"]:
            bkgnd_spec = bkgnd_spec[1]
        bkgnd_spec = read_instr_spectrum(bkgnd_spec[0], bkgnd_spec[1])
        out_spec += generate_channel_spectrum(
            bkgnd_spec, exp_time, bkgnd_area, prng=prng
        )
    if ptsrc_bkgnd:
        mylog.info("Adding in background from unresolved point-sources.")
        bkgnd_nH = float(soxs_cfg.get("soxs", "bkgnd_nH"))
        absorb_model = soxs_cfg.get("soxs", "bkgnd_absorb_model")
        spec_plaw = BackgroundSpectrum.from_powerlaw(
            1.52, 0.0, 2.0e-7, emin=0.01, emax=10.0, nbins=300000
        )
        spec_plaw.apply_foreground_absorption(bkgnd_nH, model=absorb_model)
        cspec_plaw = ConvolvedSpectrum.convolve(spec_plaw.to_spectrum(fov), arf)
        out_spec += rmf.convolve_spectrum(cspec_plaw, exp_time, prng=prng)

    bins = (np.arange(rmf.n_ch) + rmf.cmin).astype("int32")

    _write_spectrum(
        bins,
        out_spec,
        exp_time,
        rmf.header["CHANTYPE"],
        event_params,
        out_file,
        overwrite=overwrite,
    )


def simple_event_list(
    input_events,
    out_file,
    exp_time,
    instrument,
    overwrite=False,
    use_gal_coords=False,
    prng=None,
):
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.time import Time, TimeDelta
    from pathlib import PurePath

    if not out_file.endswith(".fits"):
        out_file += ".fits"
    mylog.info("Making simple observation of source in %s.", out_file)
    exp_time = parse_value(exp_time, "s")
    prng = parse_prng(prng)
    source_list, parameters = make_source_list(input_events)

    try:
        instrument_spec = instrument_registry[instrument]
    except KeyError:
        raise KeyError(f"Instrument {instrument} is not in the instrument registry!")
    if not instrument_spec["imaging"]:
        raise RuntimeError(
            f"Instrument '{instrument_spec['name']}' is not "
            f"designed for imaging observations!"
        )

    arf_file = get_data_file(instrument_spec["arf"])
    rmf_file = get_data_file(instrument_spec["rmf"])
    arf = AuxiliaryResponseFile(arf_file)
    rmf = RedistributionMatrixFile(rmf_file)

    event_params = {
        "exposure_time": exp_time,
        "arf": arf.filename,
        "rmf": rmf.filename,
        "channel_type": rmf.chan_type.upper(),
        "telescope": rmf.header["TELESCOP"],
        "instrument": instrument_spec["name"],
        "mission": rmf.header.get("MISSION", ""),
        "nchan": rmf.n_ch,
        "chan_lim": [rmf.cmin, rmf.cmax],
    }

    all_events = defaultdict(list)

    for i, src in enumerate(source_list):

        mylog.info("Detecting events from source %s.", parameters["src_names"][i])

        mylog.info(
            "Applying energy-dependent effective area from %s.",
            os.path.split(arf.filename)[-1],
        )
        refband = [parameters["emin"][i], parameters["emax"][i]]
        if src.src_type == "phlist":
            events = arf.detect_events_phlist(
                src.events.copy(), exp_time, parameters["flux"][i], refband, prng=prng
            )
        elif src.src_type.endswith("spectrum"):
            events = arf.detect_events_spec(src, exp_time, refband, prng=prng)

        n_evt = events["energy"].size

        if n_evt == 0:
            mylog.warning("No events were observed for this source!!!")
        else:
            # Add times to events
            events["time"] = prng.uniform(size=n_evt, low=0.0, high=exp_time)

        if n_evt > 0:
            for key in events:
                all_events[key] = np.concatenate([all_events[key], events[key]])

    if len(all_events["energy"]) == 0:
        mylog.warning(
            "No events from any of the sources in " "the catalog were detected!"
        )
        for key in ["ra", "dec", "time", event_params["channel_type"]]:
            all_events[key] = np.array([])
    else:
        # Step 4: Scatter energies with RMF
        mylog.info("Scattering energies with RMF %s.", os.path.split(rmf.filename)[-1])
        all_events = rmf.scatter_energies(all_events, prng=prng)

    mylog.info("Writing events to file %s.", out_file)

    t_begin = Time.now()
    dt = TimeDelta(event_params["exposure_time"], format="sec")
    t_end = t_begin + dt

    if use_gal_coords:
        names = ["GLON", "GLAT"]
        c = SkyCoord(all_events["ra"], all_events["dec"], unit="deg")
        lon = c.galactic.l
        lat = c.galactic.b
    else:
        names = ["RA", "DEC"]
        lon = all_events["ra"]
        lat = all_events["dec"]

    col_ra = fits.Column(name=names[0], format="E", unit="deg", array=lon)
    col_dec = fits.Column(name=names[1], format="E", unit="deg", array=lat)
    col_e = fits.Column(
        name="ENERGY", format="E", unit="eV", array=all_events["energy"] * 1000.0
    )

    chantype = event_params["channel_type"].upper()
    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"
    col_ch = fits.Column(
        name=chantype, format="1J", unit=cunit, array=all_events[chantype]
    )

    col_t = fits.Column(name="TIME", format="1D", unit="s", array=all_events["time"])

    cols = [col_e, col_ra, col_dec, col_ch, col_t]

    coldefs = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "EVENTS"

    tbhdu.header["TLMIN4"] = event_params["chan_lim"][0]
    tbhdu.header["TLMAX4"] = event_params["chan_lim"][1]
    tbhdu.header["EXPOSURE"] = event_params["exposure_time"]
    tbhdu.header["TSTART"] = 0.0
    tbhdu.header["TSTOP"] = event_params["exposure_time"]
    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["RADECSYS"] = "FK5"
    tbhdu.header["EQUINOX"] = 2000.0
    tbhdu.header["HDUCLASS"] = "OGIP"
    tbhdu.header["HDUCLAS1"] = "EVENTS"
    tbhdu.header["HDUCLAS2"] = "ACCEPTED"
    tbhdu.header["DATE"] = t_begin.tt.isot
    tbhdu.header["DATE-OBS"] = t_begin.tt.isot
    tbhdu.header["DATE-END"] = t_end.tt.isot
    tbhdu.header["RESPFILE"] = PurePath(event_params["rmf"]).parts[-1]
    tbhdu.header["PHA_BINS"] = event_params["nchan"]
    tbhdu.header["ANCRFILE"] = PurePath(event_params["arf"]).parts[-1]
    tbhdu.header["CHANTYPE"] = event_params["channel_type"]
    tbhdu.header["MISSION"] = event_params["mission"]
    tbhdu.header["TELESCOP"] = event_params["telescope"]
    tbhdu.header["INSTRUME"] = event_params["instrument"]

    start = fits.Column(name="START", format="1D", unit="s", array=np.array([0.0]))
    stop = fits.Column(
        name="STOP",
        format="1D",
        unit="s",
        array=np.array([event_params["exposure_time"]]),
    )

    tbhdu_gti = fits.BinTableHDU.from_columns([start, stop])
    tbhdu_gti.name = "STDGTI"
    tbhdu_gti.header["TSTART"] = 0.0
    tbhdu_gti.header["TSTOP"] = event_params["exposure_time"]
    tbhdu_gti.header["HDUCLASS"] = "OGIP"
    tbhdu_gti.header["HDUCLAS1"] = "GTI"
    tbhdu_gti.header["HDUCLAS2"] = "STANDARD"
    tbhdu_gti.header["RADECSYS"] = "FK5"
    tbhdu_gti.header["EQUINOX"] = 2000.0
    tbhdu_gti.header["DATE"] = t_begin.tt.isot
    tbhdu_gti.header["DATE-OBS"] = t_begin.tt.isot
    tbhdu_gti.header["DATE-END"] = t_end.tt.isot

    hdulist = [fits.PrimaryHDU(), tbhdu, tbhdu_gti]

    fits.HDUList(hdulist).writeto(out_file, overwrite=overwrite)
