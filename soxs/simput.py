import astropy.io.fits as pyfits
import numpy as np
from soxs.utils import parse_prng, parse_value, \
    ensure_list, ensure_numpy_array
from soxs.spatial import construct_wcs
from astropy.units import Quantity
from collections import defaultdict


def _parse_catalog_entry(src):
    src_file, src_id = src.split("[")
    src_id = src_id.strip("]")
    extname, extver = src_id.lower().strip("]").split(",")
    extver = int(extver)
    return src_file, extname, extver


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
        self.src_names = ensure_numpy_array(src_names)
        self.fluxes = ensure_numpy_array(fluxes)
        self.emin = ensure_numpy_array(emin)
        self.emax = ensure_numpy_array(emax)

    @classmethod
    def from_sources(cls, src_list):
        sc = cls([], [], [], [], [])
        src_list = ensure_list(src_list)
        for src in src_list:
            sc.append(src)
        return sc

    @classmethod
    def from_models(cls, src_names, src_models, t_exp, area, prng=None):
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
        src_names = ensure_list(src_names)
        if isinstance(src_models, (tuple, list)) and len(src_models) == 2:
            src_models = [src_models]
        if len(src_names) != len(src_models):
            raise RuntimeError(f"Number of src_names ({len(src_names)}) != "
                               f"number of src_models ({len(src_models)})!")
        sources = []
        for i, src_name in enumerate(src_names):
            src = SimputPhotonList.from_models(src_name, src_models[i][0], 
                                               src_models[i][1], t_exp, 
                                               area, prng=prng)
            sources.append(src)
        return cls(sources, src_names, 
                   [src.flux for src in sources],
                   [src.emin for src in sources],
                   [src.emax for src in sources])

    @classmethod
    def from_file(cls, filename):
        """
        Generate a SIMPUT catalog object by reading it in from
        disk. 

        Parameters
        ----------
        filename : string
            The name of the SIMPUT catalog file to read the 
            catalog and photon lists from. 
        """
        from .spectra import Spectrum
        sources = []
        f_simput = pyfits.open(filename)
        fluxes = f_simput["src_cat"].data["flux"]
        src_names = f_simput["src_cat"].data["src_name"]
        e_min = f_simput["src_cat"].data["e_min"]
        e_max = f_simput["src_cat"].data["e_max"]
        ra0 = f_simput["src_cat"].data["ra"]
        dec0 = f_simput["src_cat"].data["dec"]
        spectra = f_simput["src_cat"].data["spectrum"]
        if "IMAGE" in f_simput["src_cat"].columns.names:
            images = f_simput["src_cat"].data["image"]
        else:
            images = ["NULL"]*len(spectra)
        f_simput.close()
        for i, src in enumerate(spectra):
            src_file, extname, extver = _parse_catalog_entry(src)
            # If no file is specified, assume the catalog and
            # source are in the same file
            if src_file == "":
                fn_src = filename
            else:
                fn_src = src_file
            data = pyfits.getdata(fn_src, (extname, extver))
            if extname == "phlist":
                ra = Quantity(data["ra"], "deg")
                dec = Quantity(data["dec"], "deg")
                energy = Quantity(data["energy"], "keV")
                src = SimputPhotonList(ra, dec, energy,
                                       fluxes[i], name=src_names[i])
            elif extname == "spectrum":
                emid = data["energy"]
                de = np.diff(emid)[0]
                ebins = np.append(emid-0.5*de, emid[-1]+0.5*de)
                spec = Spectrum(ebins, data["fluxdensity"])
                if images[i].upper() != "NULL":
                    src_file, extname, extver = _parse_catalog_entry(images[i])
                    # If no file is specified, assume the catalog and
                    # source are in the same file
                    if src_file == "":
                        fn_img = filename
                    else:
                        fn_img = src_file
                    imdata, imheader = pyfits.getdata(
                        fn_img, (extname, extver), header=True)
                    imhdu = pyfits.ImageHDU(data=imdata,
                                            header=imheader)
                else:
                    imhdu = None
                src = SimputSpectrum(spec, ra0[i], dec0[i],
                                     name=src_names[i], imhdu=imhdu)
            sources.append(src)

        return cls(sources, src_names, fluxes, e_min, e_max)

    def write_catalog(self, filename, overwrite=False, src_filenames=None):
        """
        Write the SIMPUT catalog and associated sources to disk.

        Parameters
        ----------
        overwrite : boolean, optional
            Whether or not to overwrite an existing file with 
            the same name. Default: False
        """
        if src_filenames is None:
            src_filenames = {}
        src_id = []
        ra = []
        dec = []
        e_min = []
        e_max = []
        flx = []
        spectrum = []
        image = []
        timing = []
        src_name = []

        extver = defaultdict(int)

        for i, src in enumerate(self.sources):
            if src.name is None:
                name = f"soxs_src_{i}"
            else:
                name = src.name
            src_id.append(i)
            ra.append(src.ra)
            dec.append(src.dec)
            e_min.append(src.emin)
            e_max.append(src.emax)
            flx.append(src.flux)
            src_name.append(name)
            if i in src_filenames:
                fn = src_filenames[i]
            else:
                fn = filename
            extver[fn, src.src_type] += 1
            spec_fn = fn if fn != filename else ""
            spec = f"{spec_fn}[{src.src_type.upper()},{extver[fn, src.src_type]}]"
            spectrum.append(spec)
            if getattr(src, "imhdu", None):
                img_fn = fn if fn != filename else ""
                img = f"{img_fn}[IMAGE,{extver[fn, src.src_type]}]"
            else:
                img = "NULL"
            image.append(img)
            timing.append("NULL")

        src_id = np.array(src_id).astype("int32")

        col1 = pyfits.Column(name='SRC_ID', format='J', array=src_id)
        col2 = pyfits.Column(name='RA', format='D', array=ra)
        col3 = pyfits.Column(name='DEC', format='D', array=dec)
        col4 = pyfits.Column(name='E_MIN', format='D', array=e_min)
        col5 = pyfits.Column(name='E_MAX', format='D', array=e_max)
        col6 = pyfits.Column(name='FLUX', format='D', array=flx)
        col7 = pyfits.Column(name='SPECTRUM', format='80A', array=spectrum)
        col8 = pyfits.Column(name='IMAGE', format='80A', array=image)
        col9 = pyfits.Column(name='TIMING', format='80A', array=timing)
        col10 = pyfits.Column(name='SRC_NAME', format='80A', array=src_name)

        coldefs = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])

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

        wrhdu.writeto(filename, overwrite=overwrite)

        for i, src in enumerate(self.sources):
            if i in src_filenames:
                src_fn = src_filenames[i]
                append = False
            else:
                src_fn = filename
                append = True
            src.write_source(src_fn, append=append, overwrite=overwrite)

    def append(self, source):
        """
        Add a source to this catalog.

        Parameters
        ----------
        source : :class:`~soxs.simput.SimputSource`
            The source to append to this catalog.
        """
        self.sources.append(source)
        self.src_names = np.append(self.src_names, source.name)
        self.fluxes = np.append(self.fluxes, source.flux)
        self.emin = np.append(self.emin, source.emin)
        self.emax = np.append(self.emax, source.emax)


class SimputSource:
    src_type = "null"

    def __init__(self, emin, emax, flux, ra, dec, name=None):
        self.emin = emin
        self.emax = emax
        self.name = name
        self.flux = flux
        self.ra = ra
        self.dec = dec

    def _write_source(self, coldefs, filename, append, overwrite, 
                      header, imhdu=None):
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = self.src_type.upper()

        hduclas1 = "PHOTONS" if self.src_type == "phlist" else self.src_type.upper()
        tbhdu.header["HDUCLASS"] = "HEASARC/SIMPUT"
        tbhdu.header["HDUCLAS1"] = hduclas1
        tbhdu.header["HDUVERS"] = "1.1.0"
        tbhdu.header.update(header)
        if append:
            extver = 1
            f = pyfits.open(filename, mode='append')
            for hdu in f:
                if hdu.name == self.src_type.upper():
                    extver = hdu.header["EXTVER"] + 1
            tbhdu.header["EXTVER"] = extver
            f.append(tbhdu)
            if imhdu is not None:
                for hdu in f:
                    if hdu.name == "IMAGE":
                        extver = hdu.header["EXTVER"] + 1
                imhdu.header["EXTVER"] = extver
                f.append(imhdu)
            f.flush()
            f.close()
        else:
            tbhdu.header["EXTVER"] = 1
            f = [pyfits.PrimaryHDU(), tbhdu]
            if imhdu is not None:
                imhdu.header["EXTVER"] = 1
                f.append(imhdu)
            pyfits.HDUList(f).writeto(filename, overwrite=overwrite)


class SimputSpectrum(SimputSource):
    src_type = "spectrum"

    def __init__(self, spec, ra, dec, name=None, imhdu=None):
        super(SimputSpectrum, self).__init__(spec.ebins.value.min(),
                                             spec.ebins.value.max(),
                                             spec.total_flux.value, ra,
                                             dec, name=name)
        self.spec = spec
        self.imhdu = imhdu

    def write_source(self, filename, append=False, overwrite=False):
        col1 = pyfits.Column(name='ENERGY', format='E',
                             array=self.spec.emid.value)
        col2 = pyfits.Column(name='FLUXDENSITY', format='D',
                             array=self.spec.flux.value)
        cols = [col1, col2]

        coldefs = pyfits.ColDefs(cols)

        header = {"REFRA": self.ra,
                  "REFDEC": self.dec,
                  "TUNIT1": "keV",
                  "TUNIT2": "photon/(cm**2*s*keV)"}

        self._write_source(coldefs, filename, append, overwrite, header,
                           imhdu=self.imhdu)


    @classmethod
    def from_models(cls, name, spectral_model, spatial_model,
                    width, nx):
        imhdu = spatial_model.generate_image(width, nx)
        return cls(spectral_model, spatial_model.ra0,
                   spatial_model.dec0, name=name, imhdu=imhdu)


class SimputPhotonList(SimputSource):
    src_type = "phlist"

    def __init__(self, ra, dec, energy, flux, name=None):
        emin = np.asarray(energy).min()
        emax = np.asarray(energy).max()
        super(SimputPhotonList, self).__init__(
            emin, emax, flux, 0.0, 0.0, name=name)
        self.events = {"ra": ra, "dec": dec, "energy": energy}
        self.num_events = energy.size

    def __getitem__(self, item):
        return self.events[item]

    def __contains__(self, item):
        return item in self.events

    def __iter__(self):
        for key in self.events:
            yield key

    @classmethod
    def from_models(cls, name, spectral_model, spatial_model,
                    t_exp, area, prng=None):
        """
        Generate a single photon list from a spectral and a spatial
        model. 

        Parameters
        ----------
        name : string
            The name of the photon list. 
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
        prng = parse_prng(prng)
        t_exp = parse_value(t_exp, "s")
        area = parse_value(area, "cm**2")
        e = spectral_model.generate_energies(t_exp, area, prng=prng)
        ra, dec = spatial_model.generate_coords(e.size, prng=prng)
        return cls(ra, dec, e, e.flux.value, name=name)

    def write_source(self, filename, append=False, overwrite=False):
        """
        Write the photon list to disk, attaching it to a SIMPUT catalog.

        Parameters
        ----------
        filename : string
            The filename to write the photon list to. If it does not exist, 
            it will be created. Otherwise, the behavior depends on the 
            values of the *append* and *overwrite* keyword arguments.
        append : boolean, optional
            If True, append to an existing SIMPUT catalog. If False,
            overwrite an existing SIMPUT catalog with this name (if it
            exists) and create a new one with this photon list. Default: False
        overwrite : boolean, optional
            Whether or not to overwrite an existing file with 
            the same name. Default: False
        """
        col1 = pyfits.Column(name='ENERGY', format='E', array=self["energy"].value)
        col2 = pyfits.Column(name='RA', format='D', array=self["ra"].value)
        col3 = pyfits.Column(name='DEC', format='D', array=self["dec"].value)
        cols = [col1, col2, col3]

        coldefs = pyfits.ColDefs(cols)

        header = {"REFRA": 0.0,
                  "REFDEC": 0.0,
                  "TUNIT1": "keV",
                  "TUNIT2": "deg",
                  "TUNIT3": "deg"}

        self._write_source(coldefs, filename, append, overwrite, header)

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
