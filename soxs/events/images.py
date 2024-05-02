import numpy as np
from astropy import wcs
from astropy.io import fits
from tqdm.auto import tqdm

from soxs.instrument_registry import instrument_registry
from soxs.utils import create_region, get_rot_mat, mylog, parse_value

coord_types = {"sky": ("X", "Y", 2, 3), "det": ("DETX", "DETY", 6, 7)}


def make_image(
    evt_file,
    coord_type="sky",
    emin=None,
    emax=None,
    tmin=None,
    tmax=None,
    bands=None,
    expmap_file=None,
    reblock=1,
    width=None,
):
    r"""
    Generate an image by binning X-ray counts.

    Parameters
    ----------
    evt_file : string or :class:`~astropy.io.fits.HDUList`
        The name of the input event file to read, or an HDUList
        object.
    coord_type : string, optional
        The type of coordinate to bin into an image.
        Can be "sky" or "det". Default: "sky"
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the photons to put in the image, in keV.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the photons to put in the image, in keV.
    tmin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the events to be included, in seconds.
        Default is the earliest time available.
    tmax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the events to be included, in seconds.
        Default is the latest time available.
    bands : list of tuples, optional
        A list of energy bands to restrict the counts used to make the
        image, in the form of [(emin1, emax1), (emin2, emax2), ...].
        Used as an alternative to emin and emax. Default: None
    expmap_file : string, optional
        Supply an exposure map file to divide this image by
        to get a flux map. Default: None
    reblock : integer, optional
        Change this value to reblock the image to larger
        or small pixel sizes. Only supported for
        sky coordinates. Default: 1
    """
    if bands is not None:
        bands = [
            (parse_value(b[0], "keV") * 1000.0, parse_value(b[1], "keV") * 1000.0)
            for b in bands
        ]
    else:
        if emin is None:
            emin = 0.0
        else:
            emin = parse_value(emin, "keV")
        emin *= 1000.0
        if emax is None:
            emax = 100.0
        else:
            emax = parse_value(emax, "keV")
        emax *= 1000.0
    if coord_type == "det" and reblock != 1:
        raise RuntimeError(
            "Reblocking images is not supported for detector coordinates!"
        )
    if isinstance(evt_file, fits.HDUList):
        ehdu = evt_file["EVENTS"]
    else:
        ehdu = fits.open(evt_file)["EVENTS"]
    e = ehdu.data["ENERGY"]
    t = ehdu.data["TIME"]
    if tmin is None:
        tmin = 0.0
    else:
        tmin = parse_value(tmin, "s")
    if tmax is None:
        tmax = ehdu.header["EXPOSURE"]
    else:
        tmax = parse_value(tmax, "s")
    if bands is not None:
        idxs = False
        for band in bands:
            idxs |= np.logical_and(e > band[0], e < band[1])
    else:
        idxs = np.logical_and(e > emin, e < emax)
    idxs &= np.logical_and(t > tmin, t < tmax)
    xcoord, ycoord, xcol, ycol = coord_types[coord_type]
    x = ehdu.data[xcoord][idxs]
    y = ehdu.data[ycoord][idxs]
    xmin = ehdu.header[f"TLMIN{xcol}"]
    ymin = ehdu.header[f"TLMIN{ycol}"]
    xmax = ehdu.header[f"TLMAX{xcol}"]
    ymax = ehdu.header[f"TLMAX{ycol}"]
    if coord_type == "sky":
        xctr = ehdu.header[f"TCRVL{xcol}"]
        yctr = ehdu.header[f"TCRVL{ycol}"]
        if width is not None:
            xmin = 0.5 * (xmin + xmax) - 0.5 * width * (xmax - xmin)
            xmax = 0.5 * (xmin + xmax) + 0.5 * width * (xmax - xmin)
            ymin = 0.5 * (ymin + ymax) - 0.5 * width * (ymax - ymin)
            ymax = 0.5 * (ymin + ymax) + 0.5 * width * (ymax - ymin)
        xdel = ehdu.header[f"TCDLT{xcol}"] * reblock
        ydel = ehdu.header[f"TCDLT{ycol}"] * reblock

    nx = int(int(xmax - xmin) // reblock)
    ny = int(int(ymax - ymin) // reblock)

    xbins = np.linspace(xmin, xmax, nx + 1, endpoint=True)
    ybins = np.linspace(ymin, ymax, ny + 1, endpoint=True)

    H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])

    if expmap_file is not None:
        if coord_type == "det":
            raise RuntimeError(
                "Cannot divide by an exposure map for images "
                "binned in detector coordinates!"
            )
        with fits.open(expmap_file) as f:
            if f["EXPMAP"].shape != (nx, ny):
                raise RuntimeError(
                    "Exposure map and image do not have the same shape!!"
                )
            with np.errstate(invalid="ignore", divide="ignore"):
                H /= f["EXPMAP"].data.T
            H[np.isinf(H)] = 0.0
            H = np.nan_to_num(H)
            H[H < 0.0] = 0.0

    hdu = fits.PrimaryHDU(H.T)

    if coord_type == "sky":
        hdu.header["MTYPE1"] = "EQPOS"
        hdu.header["MFORM1"] = "RA,DEC"
        hdu.header["CTYPE1"] = "RA---TAN"
        hdu.header["CTYPE2"] = "DEC--TAN"
        hdu.header["CRVAL1"] = xctr
        hdu.header["CRVAL2"] = yctr
        hdu.header["CUNIT1"] = "deg"
        hdu.header["CUNIT2"] = "deg"
        hdu.header["CDELT1"] = xdel
        hdu.header["CDELT2"] = ydel
        hdu.header["CRPIX1"] = 0.5 * (nx + 1)
        hdu.header["CRPIX2"] = 0.5 * (ny + 1)
    else:
        hdu.header["CUNIT1"] = "pixel"
        hdu.header["CUNIT2"] = "pixel"

    hdu.header["ENERGYLO"] = emin
    hdu.header["ENERGYHI"] = emax
    hdu.header["EXPOSURE"] = tmax - tmin
    hdu.name = "IMAGE"

    return hdu


def write_image(
    evt_file,
    out_file,
    coord_type="sky",
    emin=None,
    emax=None,
    tmin=None,
    tmax=None,
    bands=None,
    overwrite=False,
    expmap_file=None,
    reblock=1,
):
    r"""
    Generate a image by binning X-ray counts and write
    it to a FITS file.

    Parameters
    ----------
    evt_file : string
        The name of the input event file to read.
    out_file : string
        The name of the image file to write.
    coord_type : string, optional
        The type of coordinate to bin into an image.
        Can be "sky" or "det". Default: "sky"
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the photons to put in the image, in keV.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the photons to put in the image, in keV.
    tmin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the events to be included, in seconds.
        Default is the earliest time available.
    tmax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the events to be included, in seconds.
        Default is the latest time available.
    bands : list of tuples, optional
        A list of energy bands to restrict the counts used to make the
        image, in the form of [(emin1, emax1), (emin2, emax2), ...].
        Used as an alternative to emin and emax. Default: None
    overwrite : boolean, optional
        Whether to overwrite an existing file with
        the same name. Default: False
    expmap_file : string, optional
        Supply an exposure map file to divide this image by
        to get a flux map. Default: None
    reblock : integer, optional
        Change this value to reblock the image to larger
        or small pixel sizes. Only supported for
        sky coordinates. Default: 1
    """
    hdu = make_image(
        evt_file,
        coord_type=coord_type,
        emin=emin,
        emax=emax,
        tmin=tmin,
        tmax=tmax,
        bands=bands,
        expmap_file=expmap_file,
        reblock=reblock,
    )
    hdu.writeto(out_file, overwrite=overwrite)


def make_exposure_map(
    event_file,
    expmap_file,
    energy,
    weights=None,
    asol_file=None,
    normalize=True,
    overwrite=False,
    reblock=1,
    nhistx=16,
    nhisty=16,
):
    """
    Make an exposure map for a SOXS event file, and optionally write
    an aspect solution file. The exposure map will be created by
    binning an aspect histogram over the range of the aspect solution.

    Parameters
    ----------
    event_file : string
        The path to the event file to use for making the exposure map.
    expmap_file : string
        The path to write the exposure map file to.
    energy : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, or NumPy array
        The energy in keV to use when computing the exposure map, or
        a set of energies to be used with the *weights* parameter. If
        providing a set, it must be in keV.
    weights : array-like, optional
        The weights to use with a set of energies given in the
        *energy* parameter. Used to create a more accurate exposure
        map weighted by a range of energies. Default: None
    asol_file : string, optional
        The path to write the aspect solution file to, if desired.
        Default: None
    normalize : boolean, optional
        If True, the exposure map will be divided by the exposure time
        so that the map's units are cm**2. Default: True
    overwrite : boolean, optional
        Whether to overwrite an existing file. Default: False
    reblock : integer, optional
        Supply an integer power of 2 here to make an exposure map
        with a different binning. Default: 1
    nhistx : integer, optional
        The number of bins in the aspect histogram in the DETX
        direction. Default: 16
    nhisty : integer, optional
        The number of bins in the aspect histogram in the DETY
        direction. Default: 16
    order : integer, optional
        The interpolation order to use when making the exposure map.
        Default: 1
    """
    from scipy.ndimage import rotate

    from soxs.instrument import perform_dither
    from soxs.response import AuxiliaryResponseFile

    if isinstance(energy, np.ndarray) and weights is None:
        raise RuntimeError(
            "Must supply a single value for the energy if " "you do not supply weights!"
        )
    if not isinstance(energy, np.ndarray):
        energy = parse_value(energy, "keV")
    f_evt = fits.open(event_file)
    hdu = f_evt["EVENTS"]
    arf = AuxiliaryResponseFile(hdu.header["ANCRFILE"])
    exp_time = hdu.header["EXPOSURE"]
    nx = int(hdu.header["TLMAX2"] - 0.5) // 2
    ny = int(hdu.header["TLMAX3"] - 0.5) // 2
    ra0 = hdu.header["TCRVL2"]
    dec0 = hdu.header["TCRVL3"]
    xdel = hdu.header["TCDLT2"]
    ydel = hdu.header["TCDLT3"]
    x0 = hdu.header["TCRPX2"]
    y0 = hdu.header["TCRPX3"]
    xaim = hdu.header.get("AIMPT_X", 0.0)
    yaim = hdu.header.get("AIMPT_Y", 0.0)
    xaim += hdu.header.get("AIMPT_DX", 0.0)
    yaim += hdu.header.get("AIMPT_DY", 0.0)
    roll = hdu.header["ROLL_PNT"]
    instr = instrument_registry[hdu.header["INSTRUME"]]
    dither_params = {}
    if "DITHXAMP" in hdu.header:
        dither_params["x_amp"] = hdu.header["DITHXAMP"]
        dither_params["y_amp"] = hdu.header["DITHYAMP"]
        dither_params["x_period"] = hdu.header["DITHXPER"]
        dither_params["y_period"] = hdu.header["DITHYPER"]
        dither_params["plate_scale"] = ydel * 3600.0
        dither_params["dither_on"] = True
    else:
        dither_params["dither_on"] = False
    f_evt.close()

    # Create time array for aspect solution
    dt = 1.0  # Seconds
    t = np.arange(0.0, exp_time + dt, dt)

    # Construct WCS
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [x0, y0]
    w.wcs.cdelt = [xdel, ydel]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ["deg"] * 2

    # Create aspect solution if we had dithering.
    # otherwise just set the offsets to zero
    if dither_params["dither_on"]:
        x_off, y_off = perform_dither(t, dither_params)
        # Make the aspect histogram
        x_amp = dither_params["x_amp"] / dither_params["plate_scale"]
        y_amp = dither_params["y_amp"] / dither_params["plate_scale"]
        x_edges = np.linspace(-x_amp, x_amp, nhistx + 1, endpoint=True)
        y_edges = np.linspace(-y_amp, y_amp, nhisty + 1, endpoint=True)
        asphist = np.histogram2d(x_off, y_off, (x_edges, y_edges))[0]
        asphist *= dt
        x_mid = 0.5 * (x_edges[1:] + x_edges[:-1]) / reblock
        y_mid = 0.5 * (y_edges[1:] + y_edges[:-1]) / reblock
    else:
        asphist = exp_time * np.ones((1, 1))
        x_mid = np.zeros(nhistx)
        y_mid = np.zeros(nhisty)

    # Determine the effective area
    eff_area = arf.interpolate_area(energy).value
    if weights is not None:
        eff_area = np.average(eff_area, weights=weights)

    rtypes = []
    args = []
    for chip in instr["chips"]:
        rtypes.append(chip[0])
        args.append(np.array(chip[1:]) / reblock)

    xdet0 = 0.5 * (2 * nx // reblock + 1)
    ydet0 = 0.5 * (2 * ny // reblock + 1)
    xaim //= reblock
    yaim //= reblock
    dx = xdet0 - xaim - 1.0
    dy = ydet0 - yaim - 1.0

    if dither_params["dither_on"]:
        niterx = nhistx
        nitery = nhisty
    else:
        niterx = 1
        nitery = 1

    expmap = np.zeros((2 * nx // reblock, 2 * ny // reblock))
    niter = niterx * nitery
    pbar = tqdm(leave=True, total=niter, desc="Creating exposure map ")
    for i in range(niterx):
        for j in range(nitery):
            chips, _ = create_region(rtypes[0], args[0], dx + x_mid[i], dy + y_mid[j])
            for rtype, arg in zip(rtypes[1:], args[1:]):
                r, _ = create_region(rtype, arg, dx + x_mid[i], dy + y_mid[j])
                chips = chips | r
            dexp = chips.to_mask().to_image(expmap.shape).astype("float64")
            expmap += dexp * asphist[i, j]
        pbar.update(nitery)
    pbar.close()

    expmap *= eff_area
    if normalize:
        expmap /= exp_time

    if roll != 0.0:
        rotate(expmap, roll, output=expmap, reshape=False)

    expmap[expmap < 0.0] = 0.0

    map_header = {
        "EXPOSURE": exp_time,
        "MTYPE1": "EQPOS",
        "MFORM1": "RA,DEC",
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "CRVAL1": ra0,
        "CRVAL2": dec0,
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CDELT1": xdel * reblock,
        "CDELT2": ydel * reblock,
        "CRPIX1": 0.5 * (2.0 * nx // reblock + 1),
        "CRPIX2": 0.5 * (2.0 * ny // reblock + 1),
    }

    map_hdu = fits.ImageHDU(expmap, header=fits.Header(map_header))
    map_hdu.name = "EXPMAP"
    map_hdu.writeto(expmap_file, overwrite=overwrite)

    if asol_file is not None:
        if dither_params["dither_on"]:
            det = np.array([x_off, y_off])

            pix = np.dot(get_rot_mat(roll).T, det)

            ra, dec = w.wcs_pix2world(pix[0, :] + x0, pix[1, :] + y0, 1)

            col_t = fits.Column(name="time", format="D", unit="s", array=t)
            col_ra = fits.Column(name="ra", format="D", unit="deg", array=ra)
            col_dec = fits.Column(name="dec", format="D", unit="deg", array=dec)

            coldefs = fits.ColDefs([col_t, col_ra, col_dec])
            tbhdu = fits.BinTableHDU.from_columns(coldefs)
            tbhdu.name = "ASPSOL"
            tbhdu.header["EXPOSURE"] = exp_time

            hdulist = [fits.PrimaryHDU(), tbhdu]

            fits.HDUList(hdulist).writeto(asol_file, overwrite=overwrite)

        else:
            mylog.warning(
                "Refusing to write an aspect solution file because "
                "there was no dithering."
            )


def plot_image(
    img_file,
    hdu="IMAGE",
    stretch="linear",
    vmin=None,
    vmax=None,
    facecolor="black",
    center=None,
    width=None,
    figsize=(10, 10),
    cmap=None,
):
    """
    Plot a FITS image created by SOXS using Matplotlib.

    Parameters
    ----------
    img_file : str
        The on-disk FITS image to plot.
    hdu : str or int, optional
        The image extension to plot. Default is "IMAGE"
    stretch : str, optional
        The stretch to apply to the colorbar scale. Options are "linear",
        "log", and "sqrt". Default: "linear"
    vmin : float, optional
        The minimum value of the colorbar. If not set, it will be the minimum
        value in the image.
    vmax : float, optional
        The maximum value of the colorbar. If not set, it will be the maximum
        value in the image.
    facecolor : str, optional
        The color of zero-valued pixels. Default: "black"
    center : array-like
        A 2-element object giving an (RA, Dec) coordinate for the center
        in degrees. If not set, the reference pixel of the image (usually
        the center) is used.
    width : float, optional
        The width of the image in degrees. If not set, the width of the
        entire image will be used.
    figsize : tuple, optional
        A 2-tuple giving the size of the image in inches, e.g. (12, 15).
        Default: (10,10)
    cmap : str, optional
        The colormap to be used. If not set, the default Matplotlib
        colormap will be used.

    Returns
    -------
    A tuple of the :class:`~matplotlib.figure.Figure` and the
    :class:`~matplotlib.axes.Axes` objects.
    """
    import matplotlib.pyplot as plt
    from astropy.visualization.wcsaxes import WCSAxes
    from astropy.wcs.utils import proj_plane_pixel_scales
    from matplotlib.colors import LogNorm, Normalize, PowerNorm

    if stretch == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    elif stretch == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif stretch == "sqrt":
        norm = PowerNorm(0.5, vmin=vmin, vmax=vmax)
    else:
        raise RuntimeError(f"'{stretch}' is not a valid stretch!")
    with fits.open(img_file) as f:
        hdu = f[hdu]
        w = wcs.WCS(hdu.header)
        pix_scale = proj_plane_pixel_scales(w)
        if center is None:
            center = w.wcs.crpix
        else:
            center = w.wcs_world2pix(center[0], center[1], 0)
        if width is None:
            dx_pix = 0.5 * hdu.shape[0]
            dy_pix = 0.5 * hdu.shape[1]
        else:
            dx_pix = width / pix_scale[0]
            dy_pix = width / pix_scale[1]
        fig = plt.figure(figsize=figsize)
        ax = WCSAxes(fig, [0.15, 0.1, 0.8, 0.8], wcs=w)
        fig.add_axes(ax)
        im = ax.imshow(hdu.data, norm=norm, cmap=cmap)
        ax.set_xlim(center[0] - 0.5 * dx_pix, center[0] + 0.5 * dx_pix)
        ax.set_ylim(center[1] - 0.5 * dy_pix, center[1] + 0.5 * dy_pix)
        ax.set_facecolor(facecolor)
        plt.colorbar(im)
    return fig, ax
