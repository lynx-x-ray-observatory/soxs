import astropy.io.fits as pyfits
import numpy as np

from xrs_tools.constants import erg_per_keV

def read_simput_phlist(simput_file):
    r"""
    Read events from a SIMPUT photon list.

    Parameters
    ----------
    simput_file : string
        The SIMPUT file to read from.

    Returns
    -------
    Two Python dictionaries:

      1. NumPy arrays of the positions and energies of the events.
      2. A set of parameters.
    """
    events = {}
    parameters = {}
    f_simput = pyfits.open(simput_file)
    parameters["flux"] = f_simput["src_cat"].data["flux"][0]
    parameters["emin"] = f_simput["src_cat"].data["e_min"][0]
    parameters["emax"] = f_simput["src_cat"].data["e_max"][0]
    phlist_file = f_simput["src_cat"].data["spectrum"][0].split("[")[0]
    f_simput.close()
    f_phlist = pyfits.open(phlist_file)
    for key in ["ra", "dec", "energy"]:
        events[key] = f_phlist["phlist"].data[key]
    f_phlist.close()
    return events, parameters

def write_simput_phlist(prefix, exp_time, area, ra, dec, energy, 
                        time=None, clobber=False, emin=None, emax=None):
    r"""
    Write events to a SIMPUT photon list.

    Parameters
    ----------
    prefix : string
        The filename prefix.
    exp_time : float
        The exposure time in seconds.
    area : float
        The effective area in cm^2.
    ra : NumPy array
        The right ascension of the photons, in degrees.
    dec : NumPy array
        The declination of the photons, in degrees.
    energy : NumPy array
        The energy of the photons, in keV.
    time : NumPy array, optional
        The arrival times of the photons, in seconds. Not included if None. 
    clobber : boolean, optional
        Set to True to overwrite previous files.
    emin : float, optional
        The minimum energy of the photons to save in keV.
    emax : float, optional
        The maximum energy of the photons to save in keV.
    """
    if emin is None:
        emin = energy.min()
    if emax is None:
        emax = energy.max()

    idxs = np.logical_and(energy >= emin, energy <= emax)
    flux = np.sum(energy[idxs])*erg_per_keV / exp_time / area

    col1 = pyfits.Column(name='ENERGY', format='E', array=energy[idxs])
    col2 = pyfits.Column(name='RA', format='D', array=ra[idxs])
    col3 = pyfits.Column(name='DEC', format='D', array=dec[idxs])
    cols = [col1, col2, col3]

    if time is not None:
        col4 = pyfits.Column(name='DEC', format='D', array=dec[idxs])
        cols.append(col4)

    coldefs = pyfits.ColDefs(cols)

    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.update_ext_name("PHLIST")

    tbhdu.header["HDUCLASS"] = "HEASARC/SIMPUT"
    tbhdu.header["HDUCLAS1"] = "PHOTONS"
    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["EXTVER"] = 1
    tbhdu.header["REFRA"] = 0.0
    tbhdu.header["REFDEC"] = 0.0
    tbhdu.header["TUNIT1"] = "keV"
    tbhdu.header["TUNIT2"] = "deg"
    tbhdu.header["TUNIT3"] = "deg"

    phfile = prefix+"_phlist.fits"

    tbhdu.writeto(phfile, clobber=clobber)

    col1 = pyfits.Column(name='SRC_ID', format='J', array=np.array([1]).astype("int32"))
    col2 = pyfits.Column(name='RA', format='D', array=np.array([0.0]))
    col3 = pyfits.Column(name='DEC', format='D', array=np.array([0.0]))
    col4 = pyfits.Column(name='E_MIN', format='D', array=np.array([float(emin)]))
    col5 = pyfits.Column(name='E_MAX', format='D', array=np.array([float(emax)]))
    col6 = pyfits.Column(name='FLUX', format='D', array=np.array([flux]))
    col7 = pyfits.Column(name='SPECTRUM', format='80A', array=np.array([phfile+"[PHLIST,1]"]))
    col8 = pyfits.Column(name='IMAGE', format='80A', array=np.array([phfile+"[PHLIST,1]"]))
    col9 = pyfits.Column(name='SRC_NAME', format='80A', array=np.array(["xrs_tools"]))

    coldefs = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])

    wrhdu = pyfits.BinTableHDU.from_columns(coldefs)
    wrhdu.update_ext_name("SRC_CAT")

    wrhdu.header["HDUCLASS"] = "HEASARC"
    wrhdu.header["HDUCLAS1"] = "SIMPUT"
    wrhdu.header["HDUCLAS2"] = "SRC_CAT"
    wrhdu.header["HDUVERS"] = "1.1.0"
    wrhdu.header["RADECSYS"] = "FK5"
    wrhdu.header["EQUINOX"] = 2000.0
    wrhdu.header["TUNIT2"] = "deg"
    wrhdu.header["TUNIT3"] = "deg"
    wrhdu.header["TUNIT4"] = "keV"
    wrhdu.header["TUNIT5"] = "keV"
    wrhdu.header["TUNIT6"] = "erg/s/cm**2"

    simputfile = prefix+"_simput.fits"

    wrhdu.writeto(simputfile, clobber=clobber)
