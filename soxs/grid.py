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


def make_grid_image(evtfile_list, out_file, emin=None, emax=None,
                    use_expmap=False, expmap_file=None,
                    expmap_energy=None, overwrite=False, reblock=1):
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
                              overwrite=overwrite, reblock=reblock)
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

