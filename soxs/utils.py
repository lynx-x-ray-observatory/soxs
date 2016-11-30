import os
import logging
import astropy.io.fits as pyfits
import numpy as np
from copy import copy

soxsLogger = logging.getLogger("soxs")

ufstring = "%(name)-3s : [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s : [%(levelname)-18s] %(asctime)s %(message)s"

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
    tbhdu.header["RESPFILE"] = os.path.split(parameters["rmf"])[-1]
    tbhdu.header["PHA_BINS"] = parameters["nchan"]
    tbhdu.header["ANCRFILE"] = os.path.split(parameters["arf"])[-1]
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

def get_wabs_absorb(e, nH):
    from soxs.spectra import wabs_cross_section
    sigma = wabs_cross_section(e)
    return np.exp(-nH*1.0e22*sigma)

def write_spectrum(evtfile, specfile, clobber=False):
    r"""
    Bin event energies into a spectrum and write it to a FITS binary table. Can bin
    on energy or channel. In the latter case, the spectral binning will be determined by
    the RMF binning.

    Parameters
    ----------
    specfile : string
        The name of the FITS file to be written.
    bin_type : string, optional
        Bin on "energy" or "channel". If an RMF is detected, channel information will be
        imported from it. 
    emin : float, optional
        The minimum energy of the spectral bins in keV. Only used if binning without an RMF.
    emax : float, optional
        The maximum energy of the spectral bins in keV. Only used if binning without an RMF.
    nchan : integer, optional
        The number of channels. Only used if binning without an RMF.
    """
    from soxs.instrument import RedistributionMatrixFile
    f = pyfits.open(evtfile)
    spectype = f["EVENTS"].header["CHANTYPE"]
    rmf = RedistributionMatrixFile(f["EVENTS"].header["RESPFILE"])
    minlength = rmf.n_ch
    if rmf.cmin == 1:
        minlength += 1
    spec = np.bincount(f["EVENTS"].data[spectype], minlength=minlength)
    if rmf.cmin == 1:
        spec = spec[1:]
    bins = (np.arange(rmf.n_ch)+rmf.cmin).astype("int32")

    col1 = pyfits.Column(name='CHANNEL', format='1J', array=bins)
    col2 = pyfits.Column(name=spectype.upper(), format='1D', array=bins.astype("float64"))
    col3 = pyfits.Column(name='COUNTS', format='1J', array=spec.astype("int32"))
    col4 = pyfits.Column(name='COUNT_RATE', format='1D', array=spec/f["EVENTS"].header["EXPOSURE"])

    coldefs = pyfits.ColDefs([col1, col2, col3, col4])

    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.update_ext_name("SPECTRUM")

    tbhdu.header["DETCHANS"] = spec.size
    tbhdu.header["TOTCTS"] = spec.sum()
    tbhdu.header["EXPOSURE"] = f["EVENTS"].header["EXPOSURE"]
    tbhdu.header["LIVETIME"] = f["EVENTS"].header["EXPOSURE"]
    tbhdu.header["CONTENT"] = spectype
    tbhdu.header["HDUCLASS"] = "OGIP"
    tbhdu.header["HDUCLAS1"] = "SPECTRUM"
    tbhdu.header["HDUCLAS2"] = "TOTAL"
    tbhdu.header["HDUCLAS3"] = "TYPE:I"
    tbhdu.header["HDUCLAS4"] = "COUNT"
    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["HDUVERS1"] = "1.1.0"
    tbhdu.header["CHANTYPE"] = spectype
    tbhdu.header["BACKFILE"] = "none"
    tbhdu.header["CORRFILE"] = "none"
    tbhdu.header["POISSERR"] = True
    for key in ["RESPFILE", "ANCRFILE", "MISSION", "TELESCOP", "INSTRUME"]:
        tbhdu.header[key] = f["EVENTS"].header[key]
    tbhdu.header["AREASCAL"] = 1.0
    tbhdu.header["CORRSCAL"] = 0.0
    tbhdu.header["BACKSCAL"] = 1.0

    f.close()

    hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

    hdulist.writeto(specfile, clobber=clobber)

def convert_rmf(rmffile):

    f = pyfits.open(rmffile)

    names = [ff.name for ff in f]
    idx = -1
    for i, name in enumerate(names):
        if "MATRIX" in name:
            idx = i
            break

    new_f = copy(f)

    matrix = new_f[names[idx]]
    fchan = matrix.data["F_CHAN"]
    nchan = matrix.data["N_CHAN"]
    m = matrix.data["MATRIX"]

    fchan_new = pyfits.Column(name = 'F_CHAN', format='PJ()',
                              array=np.array([fc[fc > 0] for fc in fchan], dtype=np.object))
    nchan_new = pyfits.Column(name = 'N_CHAN', format='PJ()',
                              array=np.array([nc[nc > 0] for nc in nchan], dtype=np.object))
    m_new = pyfits.Column(name = 'MATRIX', format='PE()',
                          array=np.array([mm[mm > 0] for mm in m], dtype=np.object))

    matrix_new = pyfits.BinTableHDU.from_columns([matrix.columns["ENERG_LO"],
                                                  matrix.columns["ENERG_HI"],
                                                  matrix.columns["N_GRP"],
                                                  fchan_new, nchan_new,
                                                  m_new])

    matrix_new.update_ext_name("SPECRESP MATRIX")

    new_f.pop(idx)

    new_f.append(matrix_new)

    header_keys = ["HDUCLASS", "HDUCLAS1", "HDUVERS1", "HDUCLAS2", "HDUVERS2", "HDUCLAS3",
                   "TELESCOP", "INSTRUME", "DETNAM", "FILTER", "DETCHANS", "CHANTYPE",
                   "HIERARCH LO_THRESH", "LO_THRES", "RMFVERSN", "TLMIN4", "TLMAX4"]

    for key in header_keys:
        if key in matrix.header:
            matrix_new.header[key] = matrix.header[key]

    new_f.writeto(os.path.split(rmffile)[-1], clobber=True)

def bin_profile(x, y, x0, y0, rmin, rmax, nbins):
    r = np.sqrt((x-x0)**2+(y-y0)**2)
    rbin = np.linspace(rmin, rmax, nbins+1)
    S, _ = np.histogram(r, bins=rbin)
    return rbin, S