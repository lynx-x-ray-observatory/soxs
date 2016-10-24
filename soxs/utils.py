import os
import logging
import astropy.io.fits as pyfits
import numpy as np
import astropy.wcs as pywcs

soxsLogger = logging.getLogger("soxs")

ufstring = "%(name)-3s: [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s: [%(levelname)-18s] %(asctime)s %(message)s"

soxs_sh = logging.StreamHandler()
# create formatter and add it to the handlers
formatter = logging.Formatter(ufstring)
soxs_sh.setFormatter(formatter)
# add the handler to the logger
soxsLogger.addHandler(soxs_sh)
soxsLogger.setLevel('INFO')
soxsLogger.propagate = False

mylog = soxsLogger

mylog.setLevel('INFO')

soxs_path = os.path.abspath(os.path.dirname(__file__))
soxs_files_path = os.path.join(soxs_path, "files")

def check_file_location(fn, subdir):
    if os.path.exists(fn):
        return os.path.abspath(fn)
    else:
        sto_fn = os.path.join(soxs_path, subdir, fn)
        if os.path.exists(sto_fn):
            return sto_fn
    raise IOError("Could not find file %s!" % fn)

def iterable(obj):
    """
    Grabbed from Python Cookbook / matploblib.cbook.  Returns true/false for
    *obj* iterable.
    """
    try: len(obj)
    except: return False
    return True

def ensure_list(obj):
    """
    This function ensures that *obj* is a list.  Typically used to convert a
    string to a list, for instance ensuring the *fields* as an argument is a
    list.
    """
    if obj is None:
        return [obj]
    if not isinstance(obj, list):
        return [obj]
    return obj

def ensure_numpy_array(obj):
    """
    This function ensures that *obj* is a numpy array. Typically used to
    convert scalar, list or tuple argument passed to functions using Cython.
    """
    if isinstance(obj, np.ndarray):
        if obj.shape == ():
            return np.array([obj])
        # We cast to ndarray to catch ndarray subclasses
        return np.array(obj)
    elif isinstance(obj, (list, tuple)):
        return np.asarray(obj)
    else:
        return np.asarray([obj])

one_arcsec = 1.0/3600.0

def construct_wcs(ra0, dec0):
    w = pywcs.WCS(naxis=2)
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [0.0]*2
    w.wcs.cdelt = [-one_arcsec, one_arcsec]
    w.wcs.ctype = ["RA---TAN","DEC--TAN"]
    w.wcs.cunit = ["deg"]*2
    return w

def write_event_file(events, parameters, filename, clobber=False):
    from astropy.time import Time, TimeDelta

    t_begin = Time.now()
    dt = TimeDelta(parameters["exposure_time"], format='sec')
    t_end = t_begin + dt

    col_x = pyfits.Column(name='X', format='D', unit='pixel', array=events["xpix"])
    col_y = pyfits.Column(name='Y', format='D', unit='pixel', array=events["ypix"])
    col_e = pyfits.Column(name='ENERGY', format='E', unit='eV', array=events["energy"]*1000.)
    col_cx = pyfits.Column(name='CHIPX', format='D', unit='pixel', array=events["chipx"])
    col_cy = pyfits.Column(name='CHIPY', format='D', unit='pixel', array=events["chipy"])
    col_dx = pyfits.Column(name='DETX', format='D', unit='pixel', array=events["detx"])
    col_dy = pyfits.Column(name='DETY', format='D', unit='pixel', array=events["dety"])

    chantype = parameters["channel_type"]
    if chantype == "PHA":
        cunit = "adu"
    elif chantype == "PI":
        cunit = "Chan"
    col_ch = pyfits.Column(name=chantype.upper(), format='1J', unit=cunit, array=events[chantype])

    col_t = pyfits.Column(name="TIME", format='1D', unit='s', array=events['time'])

    cols = [col_e, col_x, col_y, col_ch, col_t, col_cx, col_cy, col_dx, col_dy]

    coldefs = pyfits.ColDefs(cols)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.update_ext_name("EVENTS")

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
    tbhdu.header["TLMIN2"] = 0.5
    tbhdu.header["TLMIN3"] = 0.5
    tbhdu.header["TLMAX2"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMAX3"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMIN4"] = parameters["chan_lim"][0]
    tbhdu.header["TLMAX4"] = parameters["chan_lim"][1]
    tbhdu.header["TLMIN6"] = 0.5
    tbhdu.header["TLMAX6"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMIN7"] = 0.5
    tbhdu.header["TLMAX7"] = parameters["num_pixels"]+0.5
    tbhdu.header["TLMIN8"] = 1.0-parameters["pix_center"][0]
    tbhdu.header["TLMAX8"] = parameters["num_pixels"]-parameters["pix_center"][0]
    tbhdu.header["TLMIN9"] = 1.0-parameters["pix_center"][1]
    tbhdu.header["TLMAX9"] = parameters["num_pixels"]-parameters["pix_center"][1]
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
    tbhdu.header["RESPFILE"] = parameters["rmf"]
    tbhdu.header["PHA_BINS"] = parameters["nchan"]
    tbhdu.header["ANCRFILE"] = parameters["arf"]
    tbhdu.header["CHANTYPE"] = parameters["channel_type"]
    tbhdu.header["MISSION"] = parameters["mission"]
    tbhdu.header["TELESCOP"] = parameters["telescope"]
    tbhdu.header["INSTRUME"] = parameters["instrument"]
    tbhdu.header["RA_PNT"] = parameters["sky_center"][0]
    tbhdu.header["DEC_PNT"] = parameters["sky_center"][1]
    tbhdu.header["ROLL_PNT"] = parameters["roll_angle"]

    start = pyfits.Column(name='START', format='1D', unit='s',
                          array=np.array([0.0]))
    stop = pyfits.Column(name='STOP', format='1D', unit='s',
                         array=np.array([parameters["exposure_time"]]))

    tbhdu_gti = pyfits.BinTableHDU.from_columns([start,stop])
    tbhdu_gti.update_ext_name("STDGTI")
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

    hdulist = [pyfits.PrimaryHDU(), tbhdu, tbhdu_gti]

    pyfits.HDUList(hdulist).writeto(filename, clobber=clobber)