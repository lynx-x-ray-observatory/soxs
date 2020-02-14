import astropy.io.fits as pyfits
import numpy as np
import os
from soxs.utils import parse_prng, parse_value, \
    ensure_list, mylog, ensure_numpy_array
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
    sc = SimputCatalog.from_file(simput_file)
    parameters = {"emin": sc.emin,
                  "emax": sc.emax,
                  "src_names": sc.src_names,
                  "flux": sc.fluxes}
    return sc.sources, parameters


def handle_simput_catalog(simputfile, phfiles, flux, emin, emax, 
                          src_names, append, overwrite):

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
    energy = np.asarray(energy)
    ra = np.asarray(ra)
    dec = np.asarray(dec)
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


def _e_from_models(prng, t_exp, area, spectral_model, spatial_model):
    prng = parse_prng(prng)
    t_exp = parse_value(t_exp, "s")
    area = parse_value(area, "cm**2")
    e = spectral_model.generate_energies(t_exp, area, prng=prng)
    ra, dec = spatial_model.generate_coords(e.size, prng=prng)
    return ra, dec, e


class SimputCatalog:

    def __init__(self, sources, src_names, fluxes, emin, emax):
        """
        Create a SIMPUT catalog from a single or multiple sources.

        Parameters
        ----------
        sources : single or list of :class:`~soxs.simput.SimputSource` instances
            The photon list(s) to create this catalog with.
        """
        self.sources = ensure_list(sources)
        self.src_names = ensure_list(src_names)
        self.fluxes = ensure_numpy_array(fluxes)
        self.emin = ensure_numpy_array(emin)
        self.emax = ensure_numpy_array(emax)

    @classmethod
    def from_models(cls, src_name, spectral_model, spatial_model,
                    t_exp, area, prng=None):
        """
        Generate a SIMPUT catalog object and a single photon list
        from a spectral and a spatial model. 

        Parameters
        ----------
        src_name : string
            The name of the source This will be the prefix of 
            the source file which is created from the SimputSource
            object which is created here.
        spectral_model : :class:`~soxs.spectra.Spectrum`
            The spectral model to use to generate the event energies.
        spatial_model : :class:`~soxs.spatial.SpatialModel`
            The spatial model to use to generate the event coordinates.
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The effective area in cm**2. If one is creating  events for a 
            SIMPUT file, a constant should be used and it must be large 
            enough so that a sufficiently large sample is drawn for the ARF.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only 
            be specified if you have a reason to generate the same 
            set of random numbers, such as for a test. Default is None, 
            which sets the seed based on the system time.

        """
        ra, dec, e = _e_from_models(prng, t_exp, area, spectral_model, 
                                    spatial_model)
        photon_list = SimputPhotonList(src_name, ra, dec, e, e.flux)
        return cls(photon_list, src_name, e.flux, e.min(), e.max())

    @classmethod
    def from_file(cls, simput_file):
        """
        Generate a SIMPUT catalog object by reading it in from
        disk. 

        Parameters
        ----------
        simput_file : string
            The name of the SIMPUT catalog file to read the 
            catalog and photon lists from. 
        """
        sources = []
        f_simput = pyfits.open(simput_file)
        fluxes = f_simput["src_cat"].data["flux"]
        src_names = f_simput["src_cat"].data["src_name"]
        e_min = f_simput["src_cat"].data["e_min"]
        e_max = f_simput["src_cat"].data["e_max"]
        spectra = f_simput["src_cat"].data["spectrum"]
        f_simput.close()
        for i, src in enumerate(spectra):
            src_file, src_id = src.split("[")
            src_id = src_id.strip("]")
            # If no file is specified, assume the catalog and source are in the same file
            if src_file == "":
                fn_src = simput_file
            else:
                fn_src = src_file
            extname, extver = src_id.lower().strip("]").split(",")
            extver = int(extver)
            if extname == "phlist":
                data = pyfits.getdata(fn_src, (extname, extver))
                ra = Quantity(data["ra"], "deg")
                dec = Quantity(data["dec"], "deg")
                energy = Quantity(data["energy"], "keV")
                src = SimputPhotonList(src_names[i], ra, dec, energy,
                                       fluxes[i])
            sources.append(src)

        return cls(sources, src_names, fluxes, e_min, e_max)

    def write_catalog(self, overwrite=False):
        """
        Write the SIMPUT catalog and associated sources to disk.

        Parameters
        ----------
        overwrite : boolean, optional
            Whether or not to overwrite an existing file with 
            the same name. Default: False
        """
        for i, src in enumerate(self.sources):
            if i == 0:
                append = False
                mylog.info("Writing SIMPUT catalog file %s_simput.fits." % self.name)
            else:
                append = True
            if src.src_type == "phlist":
                mylog.info("Writing SIMPUT photon list file %s_phlist.fits." % src.name)
                src.write_photon_list(self.name, append=append, overwrite=overwrite)

    def append(self, source):
        """
        Add a source to this catalog.

        Parameters
        ----------
        source : :class:`~soxs.simput.SimputSource`
            The source to append to this catalog.
        """
        self.sources.append(source)


class SimputSource:
    src_type = "spectrum"

    def __init__(self, name, emin, emax, flux):
        self.emin = emin
        self.emax = emax
        self.name = name
        self.flux = flux

    def __getitem__(self, item):
        return self.events[item]

    def __setitem__(self, item, value):
        self.events[item] = value

    def __contains__(self, item):
        return item in self.events

    def __iter__(self):
        for key in self.events:
            yield key

    def pop(self, item):
        return self.events.pop(item)

    @classmethod
    def from_spectrum(cls, name, spec, ra, dec):
        hdu0 = pyfits.PrimaryHDU()
        col_e = pyfits.Column(name='energy', format='D', unit='keV', array=self.emid.value)
        col_f = pyfits.Column(name='flux', format='D', unit='photon/s/cm**2/keV', array=self.flux.value)
        col_name = pyfits.Column(name='name', format='80A', array=np.array([name]))

        coldefs = pyfits.ColDefs([col_e, col_f, col_name])
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = "SPECTRUM"
        tbhdu.header["HDUCLASS"] = "HEASARC"
        tbhdu.header["HDUCLAS1"] = "SIMPUT"
        tbhdu.header["HDUCLAS2"] = "SPECTRUM"
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header["EXTVER"] = 1

        hdulist = [hdu0, tbhdu]

        pyfits.HDUList(hdulist).writeto(filename, overwrite=overwrite)


class SimputPointSource(SimputSource):
    src_type = "ptsrc_spectrum"


class SimputPhotonList(SimputSource):
    src_type = "phlist"

    def __init__(self, name, ra, dec, energy, flux):
        super(SimputPhotonList, self).__init__(name, energy.value.min(),
                                               energy.value.max(), flux)
        self.events = {"ra": ra, "dec": dec, "energy": energy}
        self.num_events = energy.size

    @classmethod
    def from_models(cls, name, spectral_model, spatial_model,
                    t_exp, area, prng=None):
        """
        Generate a single photon list from a spectral and a spatial
        model. 

        Parameters
        ----------
        name : string
            The name of the photon list. This will also be the prefix 
            of any photon list file that is written from this photon list.
        spectral_model : :class:`~soxs.spectra.Spectrum`
            The spectral model to use to generate the event energies.
        spatial_model : :class:`~soxs.spatial.SpatialModel`
            The spatial model to use to generate the event coordinates.
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The effective area in cm**2. If one is creating 
            events for a SIMPUT file, a constant should be 
            used and it must be large enough so that a 
            sufficiently large sample is drawn for the ARF.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only 
            be specified if you have a reason to generate the same 
            set of random numbers, such as for a test. Default is None, 
            which sets the seed based on the system time. 
        """
        ra, dec, e = _e_from_models(prng, t_exp, area, spectral_model,
                                    spatial_model)
        return cls(name, ra, dec, e, e.flux)

    def write_photon_list(self, simput_prefix, append=False, overwrite=False):
        """
        Write the photon list to disk, attaching it to a SIMPUT catalog.

        Parameters
        ----------
        simput_prefix : string
            The prefix of the SIMPUT catalog file to associate this
            photon list with. If it does not exist, it will be created.
        append : boolean, optional
            If True, append to an existing SIMPUT catalog. If False,
            overwrite an existing SIMPUT catalog with this name (if it
            exists) and create a new one with this photon list. Default: False
        overwrite : boolean, optional
            Whether or not to overwrite an existing file with 
            the same name. Default: False
        """
        write_photon_list(simput_prefix, self.name, self.flux,
                          self["ra"], self["dec"], self["energy"],
                          append=append, overwrite=overwrite)

    def plot(self, center, width, s=None, c=None, marker=None, stride=1,
             emin=None, emax=None, label=None, fontsize=18, fig=None, 
             ax=None, **kwargs):
        """
        Plot event coordinates from this photon list in a scatter plot, 
        optionally restricting the photon energies which are plotted
        and using only a subset of the photons. 

        Parameters
        ----------
        center : array-like
            The RA, Dec of the center of the plot in degrees.
        width : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The width of the plot in arcminutes.
        s : integer, optional
            Size of the scatter marker in points^2.
        c : string, optional
            The color of the points.
        marker : string, optional
            The marker to use for the points in the scatter plot. Default: 'o'
        stride : integer, optional
            Plot every *stride* events. Default: 1
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy of the photons to plot. Default is
            the minimum energy in the list.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy of the photons to plot. Default is
            the maximum energy in the list.
        label : string, optional
            The label of the spectrum. Default: None
        fontsize : int
            Font size for labels and axes. Default: 18
        fig : :class:`~matplotlib.figure.Figure`, optional
            A Figure instance to plot in. Default: None, one will be
            created if not provided.
        ax : :class:`~matplotlib.axes.Axes`, optional
            An Axes instance to plot in. Default: None, one will be
            created if not provided.
        """
        import matplotlib.pyplot as plt
        from astropy.visualization.wcsaxes import WCSAxes
        if fig is None:
            fig = plt.figure(figsize=(10, 10))
        if ax is None:
            wcs = construct_wcs(center[0], center[1])
            ax = WCSAxes(fig, [0.15, 0.1, 0.8, 0.8], wcs=wcs)
            fig.add_axes(ax)
        else:
            wcs = ax.wcs
        if emin is None:
            emin = self.emin
        else:
            emin = parse_value(emin, "keV")
        if emax is None:
            emax = self.emax
        else:
            emax = parse_value(emax, "keV")
        idxs = np.logical_and(self["energy"].value >= emin, self["energy"].value <= emax)
        ra = self["ra"][idxs][::stride].value
        dec = self["dec"][idxs][::stride].value
        x, y = wcs.wcs_world2pix(ra, dec, 1)
        ax.scatter(x, y, s=s, c=c, marker=marker, label=label, **kwargs)
        x0, y0 = wcs.wcs_world2pix(center[0], center[1], 1)
        width = parse_value(width, "arcmin")*60.0
        ax.set_xlim(x0-0.5*width, x0+0.5*width)
        ax.set_ylim(y0-0.5*width, y0+0.5*width)
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
        ax.tick_params(axis='both', labelsize=fontsize)
        return fig, ax
