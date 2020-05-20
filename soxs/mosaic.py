from soxs.instrument import instrument_simulator
from soxs.events import write_image, make_exposure_map
from soxs.utils import parse_value, mylog
import numpy as np
from astropy.io import fits, ascii
from astropy import wcs


def make_mosaic_events(pointing_list, input_source, out_prefix, exp_time, 
                       instrument, overwrite=False, instr_bkgnd=True,
                       foreground=True, ptsrc_bkgnd=True, bkgnd_file=None, 
                       no_dither=False, dither_params=None, subpixel_res=False, 
                       prng=None):
    """
    Observe a source from many different pointings. 

    Parameters
    ----------
    pointing_list : filename
        A two-column ASCII table of RA and Dec pointings for each mock
        observation.
    input_source : string
        The path to the SIMPUT catalog file which contains the input
        source(s).
    out_prefix : string
        The prefix for the event files which will be generated. 
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
    t = ascii.read(pointing_list, format='commented_header', guess=False,
                   header_start=0, delimiter="\t")
    out_list = []
    for i, row in enumerate(t):
        out_file = f"{out_prefix}_{i}_evt.fits"
        out_list.append(out_file)
        instrument_simulator(input_source, out_file, exp_time, instrument,
                             (row["ra"], row["dec"]), overwrite=overwrite,
                             instr_bkgnd=instr_bkgnd, foreground=foreground,
                             ptsrc_bkgnd=ptsrc_bkgnd, bkgnd_file=bkgnd_file,
                             no_dither=no_dither, dither_params=dither_params,
                             subpixel_res=subpixel_res, prng=prng)
    t["evtfile"] = out_list
    outfile = f"{out_prefix}_event_mosaic.dat"
    mylog.info(f"Writing mosaic information to {outfile}.")
    t.write(outfile, overwrite=overwrite, delimiter="\t",
            format='ascii.commented_header')
    return outfile


def make_mosaic_image(evtfile_list, img_file, emin=None, emax=None,
                      reblock=1, use_expmap=False, expmap_energy=None, 
                      expmap_weights=None, normalize=True, nhistx=16,
                      nhisty=16, overwrite=False):
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
    from reproject.mosaicking import find_optimal_celestial_wcs, \
        reproject_and_coadd
    from reproject import reproject_interp
    t = ascii.read(evtfile_list, format='commented_header',
                   guess=False, header_start=0, delimiter="\t")

    files = []
    for row in t:
        evt_file = row["evtfile"]
        img_file = evt_file.replace("evt", "img")
        if use_expmap:
            expmap_file = evt_file.replace("evt", "expmap")
            make_exposure_map(evt_file, expmap_file, energy=expmap_energy,
                              weights=expmap_weights, normalize=normalize,
                              overwrite=overwrite, reblock=reblock, nhistx=nhistx, 
                              nhisty=nhisty)
        else:
            expmap_file = None
        write_image(evt_file, img_file, emin=emin, emax=emax,
                    overwrite=overwrite, reblock=reblock)
        files.append([img_file, expmap_file])

    img_hdus = [fits.open(fns[0], memmap=True)[0] for fns in files]
    wcs_out, shape_out = find_optimal_celestial_wcs(img_hdus)

    img, footprint = reproject_and_coadd(img_hdus, wcs_out, shape_out=shape_out,
                                         reproject_function=reproject_interp)

    if use_expmap:
        emap_hdus = [fits.open(fns[1], memmap=True)[1] for fns in files]
        emap, footprint = reproject_and_coadd(emap_hdus, wcs_out, shape_out=shape_out,
                                              reproject_function=reproject_interp)

