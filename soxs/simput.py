import astropy.io.fits as pyfits
import numpy as np
import os
from soxs.utils import parse_prng, parse_value, \
    ensure_list, mylog
from soxs.spatial import construct_wcs
from astropy.units import Quantity

def read_simput_catalog(simput_file):
    r"""
    Read events from a SIMPUT catalog. This will read 
    all of the sources in the catalog.

    Parameters
    ----------
    simput_file : string
        The SIMPUT file to read from.

    Returns
    -------
    1. Lists of dicts of NumPy arrays of the positions 
       and energies of the events from the sources.
    2. NumPy arrays of the parameters of the sources.
    """
    events = []
    parameters = {}
    simput_dir = os.path.split(os.path.abspath(simput_file))[0]
    f_simput = pyfits.open(simput_file)
    parameters["flux"] = f_simput["src_cat"].data["flux"]
    parameters["emin"] = f_simput["src_cat"].data["e_min"]
    parameters["emax"] = f_simput["src_cat"].data["e_max"]
    parameters["sources"] = f_simput["src_cat"].data["src_name"]
    phlist_files = [file.split("[")[0] for file in 
                    f_simput["src_cat"].data["spectrum"]]
    f_simput.close()
    for phlist_file in phlist_files:
        f_phlist = pyfits.open(os.path.join(simput_dir, phlist_file))
        evt = {}
        for key in ["ra", "dec", "energy"]:
            evt[key] = f_phlist["phlist"].data[key]
        f_phlist.close()
        events.append(evt)
    return events, parameters

def handle_simput_catalog(simput_prefix, phfiles, flux, emin, emax, 
                          src_names, append, overwrite):

    simputfile = simput_prefix+"_simput.fits"

    num_new_sources = len(phfiles)

    if append:
        f = pyfits.open(simputfile)
        num_sources = f["SRC_CAT"].data["SRC_ID"].size
        spec_files = f["SRC_CAT"].data["SPECTRUM"][:]
        img_files = f["SRC_CAT"].data["IMAGE"][:]
        names = []
        ph_exts = []
        for i, phfile in enumerate(phfiles):
            id = num_sources + 1 + i
            ph_ext = phfile + "[PHLIST,1]"
            if ph_ext in spec_files or ph_ext in img_files:
                raise IOError("This SIMPUT catalog already has an entry for file %s!" % phfile)
            ph_exts.append(ph_ext)
            if src_names[i] is None:
                names.append("soxs_src_%d" % id)
            else:
                names.append(src_names[i])
        src_id = np.concatenate([f["SRC_CAT"].data["SRC_ID"][:], 
                                 np.arange(num_new_sources)+num_sources+1])
        spectrum = np.concatenate([spec_files, ph_exts])
        image = np.concatenate([img_files, ph_exts])
        ra = np.concatenate([f["SRC_CAT"].data["RA"][:], np.zeros(num_new_sources)])
        dec = np.concatenate([f["SRC_CAT"].data["DEC"][:], np.zeros(num_new_sources)])
        e_min = np.concatenate([f["SRC_CAT"].data["E_MIN"][:], emin])
        e_max = np.concatenate([f["SRC_CAT"].data["E_MAX"][:], emax])
        flx = np.concatenate([f["SRC_CAT"].data["FLUX"][:], flux])
        src_name = np.concatenate([f["SRC_CAT"].data["SRC_NAME"][:], names])
        f.close()
    else:
        src_id = np.arange(num_new_sources).astype("int32")+1
        ra = np.array([0.0]*num_new_sources)
        dec = np.array([0.0]*num_new_sources)
        e_min = np.array([emin]*num_new_sources)
        e_max = np.array([emax]*num_new_sources)
        flx = np.array([flux]*num_new_sources)
        spectrum = np.array([phfile + "[PHLIST,1]" for phfile in phfiles])
        image = spectrum
        names = []
        for i, src_name in enumerate(src_names):
            if src_name is None:
                names.append("soxs_src_%d" % (i+1))
            else:
                names.append(src_name)
        src_name = np.array(names)

    col1 = pyfits.Column(name='SRC_ID', format='J', array=src_id)
    col2 = pyfits.Column(name='RA', format='D', array=ra)
    col3 = pyfits.Column(name='DEC', format='D', array=dec)
    col4 = pyfits.Column(name='E_MIN', format='D', array=e_min)
    col5 = pyfits.Column(name='E_MAX', format='D', array=e_max)
    col6 = pyfits.Column(name='FLUX', format='D', array=flx)
    col7 = pyfits.Column(name='SPECTRUM', format='80A', array=spectrum)
    col8 = pyfits.Column(name='IMAGE', format='80A', array=image)
    col9 = pyfits.Column(name='SRC_NAME', format='80A', array=src_name)

    coldefs = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])

    wrhdu = pyfits.BinTableHDU.from_columns(coldefs)
    wrhdu.name = "SRC_CAT"

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

    wrhdu.writeto(simputfile, overwrite=(overwrite or append))

def write_photon_list(simput_prefix, phlist_prefix, flux, ra, dec, energy,
                      append=False, overwrite=False):
    r"""
    Write events to a new SIMPUT photon list. It can be 
    associated with a new or existing SIMPUT catalog. 

    Parameters
    ----------
    simput_prefix : string
        The filename prefix for the SIMPUT file.
    phlist_prefix : string
        The filename prefix for the photon list file.
    flux : float
        The energy flux of all the photons, in units of 
        erg/cm**2/s.
    ra : NumPy array
        The right ascension of the photons, in degrees.
    dec : NumPy array
        The declination of the photons, in degrees.
    energy : NumPy array or array-like thing
        The energy of the photons, in keV.
    append : boolean, optional
        If True, append a new source an existing SIMPUT 
        catalog. Default: False
    overwrite : boolean, optional
        Set to True to overwrite previous files. Default: False
    """
    # Make sure these are arrays
    energy = np.array(energy)
    ra = np.array(ra)
    dec = np.array(dec)
    if hasattr(flux, "value"):
        flux = flux.value

    emin = energy.min()
    emax = energy.max()

    col1 = pyfits.Column(name='ENERGY', format='E', array=energy)
    col2 = pyfits.Column(name='RA', format='D', array=ra)
    col3 = pyfits.Column(name='DEC', format='D', array=dec)
    cols = [col1, col2, col3]

    coldefs = pyfits.ColDefs(cols)

    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = "PHLIST"

    tbhdu.header["HDUCLASS"] = "HEASARC/SIMPUT"
    tbhdu.header["HDUCLAS1"] = "PHOTONS"
    tbhdu.header["HDUVERS"] = "1.1.0"
    tbhdu.header["EXTVER"] = 1
    tbhdu.header["REFRA"] = 0.0
    tbhdu.header["REFDEC"] = 0.0
    tbhdu.header["TUNIT1"] = "keV"
    tbhdu.header["TUNIT2"] = "deg"
    tbhdu.header["TUNIT3"] = "deg"

    phfile = phlist_prefix+"_phlist.fits"

    tbhdu.writeto(phfile, overwrite=overwrite)

    handle_simput_catalog(simput_prefix, [phfile], [flux], [emin], 
                          [emax], [phlist_prefix], append, overwrite)

class SimputCatalog(object):

    @classmethod
    def from_file(cls, simput_file):
        photon_lists = []
        events, parameters = read_simput_catalog(simput_file)
        for i, params in enumerate(parameters):
            ra = Quantity(events[i]["ra"], "deg")
            dec = Quantity(events[i]["dec"], "deg")
            energy = Quantity(events[i]["energy"], "keV")
            phlist = PhotonList(params["sources"][i], ra, dec, energy,
                                params[i]["flux"])
            photon_lists.append(phlist)

        return cls(photon_lists)

    def __init__(self, photon_lists):
        self.photon_lists = ensure_list(photon_lists)

    def write_catalog(self, simput_prefix, overwrite=False):
        for i, phlist in enumerate(self.photon_lists):
            if i == 0:
                append = False
                mylog.info("Writing SIMPUT catalog file %s_simput.fits." % simput_prefix)
            else:
                append = True
            mylog.info("Writing SIMPUT photon list file %s_phlist.fits." % phlist.name)
            phlist.write_photon_list(simput_prefix, phlist.name,
                                     src_name=phlist.name, append=append,
                                     overwrite=overwrite)

class PhotonList(object):

    @classmethod
    def from_models(cls, name, spectral_model, spatial_model,
                    t_exp, area, prng=None):
        prng = parse_prng(prng)
        t_exp = parse_value(t_exp, "s")
        area = parse_value(area, "cm**2")
        e = spectral_model.generate_energies(t_exp, area, prng=prng)
        ra, dec = spatial_model.generate_sample(e.size, prng=prng)
        return cls(name, ra, dec, e, e.flux)

    def __init__(self, name, ra, dec, energy, flux):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.energy = energy
        self.emin = energy.value.min()
        self.emax = energy.value.max()
        self.flux = flux
        self.num_events = energy.size

    def write_photon_list(self, simput_prefix, phlist_prefix, append=False,
                          overwrite=False):
        write_photon_list(simput_prefix, phlist_prefix, self.flux,
                          self.ra, self.dec, self.energy,
                          append=append, overwrite=overwrite)

    def plot(self, center, width, s=None, marker=None, stride=1,
             emin=None, emax=None, label=None, legend_kwargs=None,
             fontsize=18, fig=None, ax=None, **kwargs):
        import matplotlib.pyplot as plt
        try:
            from wcsaxes import WCSAxes
        except ImportError:
            raise ImportError("Using the plot functionality for PhotonList "
                              "requires the WCSAxes package to be installed.")
        if legend_kwargs is None:
            legend_kwargs = {"fontsize": fontsize}
        if fig is None:
            fig = plt.figure(figsize=(10, 10))
        if ax is None:
            wcs = construct_wcs(center[0], center[1])
            ax = WCSAxes(fig, [0.15, 0.1, 0.8, 0.8], wcs=wcs)
            fig.add_axes(ax)
        else:
            wcs = ax.wcs
        if emin is None:
            emin = self.energy.value.min()
        if emax is None:
            emax = self.energy.value.max()
        idxs = np.logical_and(self.energy.value >= emin, self.energy.value <= emax)
        ra = self.ra[idxs][::stride]
        dec = self.dec[idxs][::stride]
        x, y = wcs.wcs_world2pix(ra, dec, 1)
        ax.scatter(x, y, s=s, marker=marker, label=label, **kwargs)
        x0, y0 = wcs.wcs_world2pix(center[0], center[1], 1)
        width = parse_value(width, "arcmin")*60.0
        ax.set_xlim(x0-0.5*width, x0+0.5*width)
        ax.set_ylim(y0-0.5*width, y0+0.5*width)
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
        ax.tick_params(axis='both', labelsize=fontsize)
        num_labels = len([line for line in ax.lines
                          if not line.get_label().startswith("_line")])
        if num_labels > 0:
            ax.legend(**legend_kwargs)
        return fig, ax
