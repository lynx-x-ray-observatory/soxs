from soxs.instrument import instrument_simulator
from soxs.events import make_exposure_map, \
    wcs_from_event_file
from soxs.utils import parse_value, mylog
import numpy as np
from astropy.io import fits, ascii
from astropy import wcs
from astropy.table import Table


def observe_grid_source(grid_spec_file, exp_time, instrument,
                        overwrite=False, instr_bkgnd=True,
                        foreground=True, ptsrc_bkgnd=True,
                        bkgnd_file=None, no_dither=False,
                        dither_params=None, subpixel_res=False, prng=None):
    """
    Observe a grid of sources produced by pyXSIM, or SOXS in its cosmological
    sources mode.

    Parameters
    ----------
    grid_spec_file : filename
        The ASCII table produced by pyXSIM or SOXS containing the specification
        of the input SIMPUT photon lists and their locations on the sky.
    exp_time : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The exposure time in seconds.
    instrument : string
        The name of the instrument to use, which picks an instrument
        specification from the instrument registry.
    overwrite : boolean, optional
        Whether or not to overwrite an existing file with the same name.
        Default: False
    instr_bkgnd : boolean, optional
        Whether or not to include the instrumental/particle background. 
        Default: True
    foreground : boolean, optional
        Whether or not to include the local foreground. 
        Default: True
    ptsrc_bkgnd : boolean, optional
        Whether or not to include the point-source background. 
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
    subpixel_res: boolean, optional
        If True, event positions are not randomized within the pixels 
        within which they are detected. Default: False
    prng : :class:`~numpy.random.RandomState` object, integer, or None
        A pseudo-random number generator. Typically will only 
        be specified if you have a reason to generate the same 
        set of random numbers, such as for a test. Default is None, 
        which sets the seed based on the system time. 
    """
    t = ascii.read(grid_spec_file, format='commented_header', guess=False, 
                   header_start=0, delimiter="\t")
    simput_file = t.meta["comments"][0].split()[-1]
    out_list = []
    for row in t:
        out_file = row["phlist"].split("_phlist.fits")[0]+"_evt.fits"
        out_list.append(out_file)
        instrument_simulator(simput_file, out_file, exp_time, instrument,
                             (row["ra"], row["dec"]), overwrite=overwrite,
                             instr_bkgnd=instr_bkgnd, foreground=foreground,
                             ptsrc_bkgnd=ptsrc_bkgnd, bkgnd_file=bkgnd_file,
                             no_dither=no_dither, dither_params=dither_params,
                             subpixel_res=subpixel_res, prng=prng,
                             source_id=row["src_id"])
    et = Table([out_list], names=["evtfile"])
    et.meta["comments"] = t.meta["comments"][2:]
    outfile = "{}_event_grid.txt".format(simput_file.split("_simput.fits")[0])
    mylog.info("Writing grid information to {}.".format(outfile))
    et.write(outfile, overwrite=overwrite,
             delimiter="\t", format='ascii.commented_header')
    return outfile


def make_grid_image(evtfile_list, img_file, emin=None, emax=None,
                    reblock=1, use_expmap=False, expmap_file=None,
                    expmap_energy=None, expmap_weights=None,
                    overwrite=False):
    """
    Make a single FITS image from a grid of observations. Optionally,
    an exposure map can be computed and a flux image may be generated.

    Parameters
    ----------
    evtfile_list : filename
        The ASCII table produced by :meth:`~soxs.grid.observe_grid_source`
        containing the information about the event files and their
        locations on the sky.
    img_file : filename
        The name of the FITS image file to be written. This name will
        also be used for the exposure map and flux files if they are
        written.
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The minimum energy of the photons to put in the image, in keV.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
        The maximum energy of the photons to put in the image, in keV.
    reblock : integer, optional
        Supply an integer power of 2 here to make an exposure map 
        with a different binning. Default: 1
    use_expmap : boolean, optional
        Whether or not to use (and potentially generate) an exposure map
        and a flux map. Default: False
    expmap_file : filename, optional
        If this is supplied, an existing exposure map file will be used
        instead of generated. Default: None
    expmap_energy : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, or NumPy array
        The energy in keV to use when computing the exposure map, or 
        a set of energies to be used with the *weights* parameter. If
        providing a set, it must be in keV.
    expmap_weights : array-like, optional
        The weights to use with a set of energies given in the
        *energy* parameter. Used to create a more accurate exposure
        map weighted by a range of energies. Default: None
    overwrite : boolean, optional
        Whether or not to overwrite an existing file with the same name.
        Default: False
    """
    from scipy.interpolate import interp2d
    t = ascii.read(evtfile_list, format='commented_header',
                   guess=False, header_start=0, delimiter="\t")

    if emin is None:
        emin = 0.0
    else:
        emin = parse_value(emin, "keV")
    emin *= 1000.
    if emax is None:
        emax = 100.0
    else:
        emax = parse_value(emax, "keV")
    emax *= 1000.

    evtfiles = t["evtfile"]

    ra0, dec0 = np.array(t.meta["comments"][0].split()[-1].split(","), dtype='float64')
    numx, numy = np.array(t.meta["comments"][1].split()[-1].split(","), dtype='int64')
    f = fits.open(evtfiles[0], memmap=True)
    exp_time = f["EVENTS"].header["EXPOSURE"]
    xmin = f["EVENTS"].header["TLMIN2"]
    ymin = f["EVENTS"].header["TLMIN3"]
    xmax = f["EVENTS"].header["TLMAX2"]
    ymax = f["EVENTS"].header["TLMAX3"]
    Lx = 0.5*(xmax-xmin)
    Ly = 0.5*(ymax-ymin)
    nx = int(Lx)//reblock
    ny = int(Ly)//reblock
    nxb = nx*numx
    nyb = ny*numy
    xdel = f["EVENTS"].header["TCDLT2"]*reblock
    ydel = f["EVENTS"].header["TCDLT3"]*reblock
    f.close()

    bigw = wcs.WCS(naxis=2)
    bigw.wcs.crval = [ra0, dec0]
    bigw.wcs.crpix = [0.5*(nxb+1), 0.5*(nyb+1)]
    bigw.wcs.cdelt = [xdel, ydel]
    bigw.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    bigw.wcs.cunit = ["deg"]*2

    xbins = np.linspace(0.5, nxb+0.5, nxb+1, endpoint=True)
    ybins = np.linspace(0.5, nyb+0.5, nyb+1, endpoint=True)
    xmid = 0.5*(xbins[1:]+xbins[:-1])
    ymid = 0.5*(ybins[1:]+ybins[:-1])

    xx = []
    yy = []
    for evtfile in evtfiles:
        f = fits.open(evtfile, memmap=True)
        e = f["EVENTS"].data["ENERGY"]
        idxs = np.logical_and(e > emin, e < emax)
        w = wcs_from_event_file(f)
        x, y = w.wcs_pix2world(f["EVENTS"].data["X"][idxs],
                               f["EVENTS"].data["Y"][idxs], 1)
        x, y = bigw.wcs_world2pix(x, y, 1)
        f.close()
        xx.append(x)
        yy.append(y)

    Cbig, _, _ = np.histogram2d(np.concatenate(xx),
                                np.concatenate(yy), bins=[xbins, ybins])

    if use_expmap:
        if expmap_file is None:
            expmap_file = evtfiles[0].split("_0_0_evt.fits")[0] + "_expmap.fits"
            make_exposure_map(evtfiles[0], expmap_file, expmap_energy,
                              weights=expmap_weights, overwrite=overwrite, 
                              reblock=reblock)
        f = fits.open(expmap_file)
        E = f["EXPMAP"].data.T
        f.close()
        Ebig = np.zeros((nxb, nyb))
        nxe, nye = E.shape
        xbe = np.linspace(0.5, nxe+0.5, nxe+1, endpoint=True)
        ybe = np.linspace(0.5, nye+0.5, nye+1, endpoint=True)
        xme = 0.5*(xbe[1:]+xbe[:-1])
        yme = 0.5*(ybe[1:]+ybe[:-1])
        xme -= 0.5*(nxe+1)
        yme -= 0.5*(nye+1)
        for evtfile in evtfiles:
            f = fits.open(evtfile, memmap=True)
            w = wcs_from_event_file(f)
            f.close()
            x0, y0 = bigw.wcs_world2pix(w.wcs.crval[0], w.wcs.crval[1], 1)
            ef = interp2d(xme+x0, yme+y0, E)
            Ebig += ef(xmid, ymid)
        with np.errstate(invalid='ignore', divide='ignore'):
            Fbig = Cbig / Ebig.T
        Fbig[np.isinf(Fbig)] = 0.0
        Fbig = np.nan_to_num(Fbig)
        Fbig[Fbig < 0.0] = 0.0

    header_keys = {"MTYPE1": "EQPOS",
                   "MFORM1": "RA,DEC",
                   "CTYPE1": "RA---TAN",
                   "CTYPE2": "DEC--TAN",
                   "CRVAL1": ra0,
                   "CRVAL2": dec0,
                   "CUNIT1": "deg",
                   "CUNIT2": "deg",
                   "CDELT1": xdel,
                   "CDELT2": ydel,
                   "CRPIX1": 0.5*(nxb+1),
                   "CRPIX2": 0.5*(nyb+1)}

    hdu = fits.PrimaryHDU(Cbig.T)
    hdu.header.update(header_keys)
    hdu.header["EXPOSURE"] = exp_time

    hdu.writeto(img_file, overwrite=overwrite)

    if use_expmap:
        hdue = fits.PrimaryHDU(Ebig)
        hdue.header.update(header_keys)
        out_exp_file = img_file.split(".")[0]+"_expmap.fits"
        hdue.writeto(out_exp_file, overwrite=overwrite)

        hduf = fits.PrimaryHDU(Fbig.T)
        hduf.header.update(header_keys)
        out_flux_file = img_file.split(".")[0]+"_flux.fits"
        hduf.writeto(out_flux_file, overwrite=overwrite)

