from soxs.instrument import RedistributionMatrixFile
from soxs.spectra import wabs_cross_section
import astropy.io.fits as pyfits
import numpy as np

def get_wabs_absorb(e, nH):
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
    f = pyfits.open(evtfile)
    spectype = f["EVENTS"].header["CHANTYPE"]
    rmf = RedistributionMatrixFile(f["EVENTS"].header["RESPFILE"])
    minlength = rmf.n_ch
    if rmf.cmin == 1: minlength += 1
    spec = np.bincount(f["EVENTS"].data[spectype], minlength=minlength)
    if rmf.cmin == 1: spec = spec[1:]
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
    for key in ["RESPFILE","ANCRFILE","MISSION","TELESCOP","INSTRUME"]:
        tbhdu.header[key] = f["EVENTS"].header[key]
    tbhdu.header["AREASCAL"] = 1.0
    tbhdu.header["CORRSCAL"] = 0.0
    tbhdu.header["BACKSCAL"] = 1.0

    f.close()

    hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), tbhdu])

    hdulist.writeto(specfile, clobber=clobber)
