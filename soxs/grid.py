from soxs.instrument import instrument_simulator
from soxs.events import make_exposure_map
from soxs.utils import parse_value, mylog
import numpy as np
from astropy import wcs
from astropy.io import fits, ascii
from astropy.table import Table


def observe_grid_source(grid_spec_file, exp_time, instrument,
                        overwrite=False, instr_bkgnd=True,
                        foreground=True, ptsrc_bkgnd=True,
                        bkgnd_file=None, no_dither=False,
                        dither_params=None, subpixel_res=False, prng=None):
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
    et = Table([out_list])
    et.meta["comments"] = t.meta["comments"][2:]
    outfile = "{}_event_grid.txt".format(simput_file.split("_simput.fits")[0])
    mylog.info("Writing grid information to {}.".format(outfile))
    et.write(outfile, overwrite=overwrite,
             delimiter="\t", format='ascii.commented_header')
    return outfile


def make_grid_image(evtfile_list, out_file, emin=None, emax=None, 
                    make_expmap=False, expmap_energy=None, 
                    overwrite=False, reblock=1):
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

    evtfiles = t["evtfiles"]

    ra0, dec0 = np.array(t.meta["comments"][1].split()[-1].split(","), dtype='float64')
    numx, numy = np.array(t.meta["comments"][2].split()[-1].split(","), dtype='int64')
    f = fits.open(evtfiles[0], memmap=True)
    exp_time = f["EVENTS"].header["EXPOSURE"]
    xmin = f["EVENTS"].header["TLMIN2"]
    ymin = f["EVENTS"].header["TLMIN3"]
    xmax = f["EVENTS"].header["TLMAX2"]
    ymax = f["EVENTS"].header["TLMAX3"]
    nx = int(xmax-xmin)//reblock
    ny = int(ymax-ymin)//reblock
    nxb = nx*numx
    nyb = ny*numy
    xdel = f["EVENTS"].header["TCDLT2"]*reblock
    ydel = f["EVENTS"].header["TCDLT3"]*reblock
    f.close()

    xbins = np.linspace(xmin, xmax, nx+1, endpoint=True)
    ybins = np.linspace(ymin, ymax, ny+1, endpoint=True)

    Hbig = np.zeros((nxb, nyb))

    for evtfile in evtfiles:
        f = fits.open(evtfile, memmap=True)
        e = f["EVENTS"].data["ENERGY"]
        idxs = np.logical_and(e > emin, e < emax)
        x = f["EVENTS"].data["X"][idxs]
        y = f["EVENTS"].data["Y"][idxs]
        H, _, _ = np.histogram2d(x, y, bins=[xbins, ybins])

    hdu = fits.PrimaryHDU(Hbig.T)

    hdu.header["MTYPE1"] = "EQPOS"
    hdu.header["MFORM1"] = "RA,DEC"
    hdu.header["CTYPE1"] = "RA---TAN"
    hdu.header["CTYPE2"] = "DEC--TAN"
    hdu.header["CRVAL1"] = ra0
    hdu.header["CRVAL2"] = dec0
    hdu.header["CUNIT1"] = "deg"
    hdu.header["CUNIT2"] = "deg"
    hdu.header["CDELT1"] = xdel
    hdu.header["CDELT2"] = ydel
    hdu.header["CRPIX1"] = 0.5*(nx+1)
    hdu.header["CRPIX2"] = 0.5*(ny+1)

    hdu.header["EXPOSURE"] = exp_time

    if make_expmap:
        expmap_file = evtfiles[0].split("_evt.fits")[0] + "_expmap.fits"
        make_exposure_map(evtfiles[0], expmap_file, expmap_energy,
                          overwrite=overwrite, reblock=reblock)

    hdu.writeto(out_file, overwrite=overwrite)
