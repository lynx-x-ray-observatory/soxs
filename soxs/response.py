import astropy.units as u
import numpy as np
from astropy import wcs
from astropy.io import fits
from tqdm.auto import tqdm

from soxs.constants import erg_per_keV
from soxs.instrument_registry import instrument_registry
from soxs.utils import (
    ensure_numpy_array,
    get_data_file,
    image_pos,
    mylog,
    parse_prng,
    parse_value,
    regrid_spectrum,
)


class AuxiliaryResponseFile:
    r"""
    A class for auxiliary response files (ARFs).

    Parameters
    ----------
    filename : string
        The filename of the ARF to be read.

    Examples
    --------
    >>> arf = AuxiliaryResponseFile("xrs_mucal_3x10_3.0eV.arf")
    """

    def __init__(self, filename):
        self.filename = get_data_file(filename)
        with fits.open(self.filename) as f:
            self.elo = f["SPECRESP"].data.field("ENERG_LO").copy()
            self.ehi = f["SPECRESP"].data.field("ENERG_HI").copy()
            self.ebins = np.append(self.elo, self.ehi[-1])
            self.de = np.diff(self.ebins)
            self.emid = 0.5 * (self.elo + self.ehi)
            self.eff_area = np.nan_to_num(f["SPECRESP"].data.field("SPECRESP")).astype(
                "float64"
            )
            self.max_area = self.eff_area.max()

    @classmethod
    def from_instrument(cls, name):
        """
        Return an :class:`~soxs.instrument.AuxiliaryResponseFile`
        object from the name of an existing instrument
        specification in SOXS.

        Parameters
        ----------
        name : string
            The name of the instrument specification to use
            to obtain the ARF object from.

        Examples
        --------
        >>> arf = soxs.AuxiliaryResponseFile.from_instrument("hdxi")
        """
        instr = instrument_registry.get(name, None)
        if instr is None:
            raise KeyError(f"Instrument '{name}' not in registry!")
        return cls(instr["arf"])

    def __str__(self):
        return self.filename

    def interpolate_area(self, energy):
        """
        Interpolate the effective area to the energies
        provided  by the supplied *energy* array.
        """
        earea = np.interp(
            np.asarray(energy), self.emid, self.eff_area, left=0.0, right=0.0
        )
        return u.Quantity(earea, "cm**2")

    def detect_events_spec(self, src, exp_time, prng=None):
        prng = parse_prng(prng)
        # This assumes linear binning for now!
        de = np.diff(src.energy)[0]
        ebins = np.append(src.energy - 0.5 * de, src.energy[-1] + 0.5 * de)
        f = regrid_spectrum(self.ebins, ebins, src.fluxdensity * de)
        N = np.cumsum(f * self.eff_area)
        n_ph = prng.poisson(lam=N[-1] * exp_time)
        randvec = prng.uniform(size=n_ph)
        randvec.sort()
        cumspec = np.insert(N, 0, 0.0)
        cumspec /= cumspec[-1]
        energy = np.interp(randvec, cumspec, self.ebins)
        if getattr(src, "imhdu", None):
            x, y = image_pos(src.imhdu.data.T, energy.size, prng)
            w = wcs.WCS(header=src.imhdu.header)
            w.wcs.crval = [src.ra, src.dec]
            ra, dec = w.wcs_pix2world(x, y, 1)
        else:
            pones = np.ones_like(energy)
            ra = src.ra * pones
            dec = src.dec * pones
        mylog.debug("%d events detected from this spectrum.", energy.size)
        return {"energy": energy, "ra": ra, "dec": dec}

    def detect_events_phlist(self, events, exp_time, flux, refband, prng=None):
        """
        Use the ARF to determine a subset of photons which
        will be detected.

        Parameters
        ----------
        events : dict of np.ndarrays
            The energies and positions of the photons.
        exp_time : float
            The exposure time in seconds.
        flux : float
            The total flux of the photons in erg/s/cm^2.
        refband : array_like
            A two-element array or list containing the limits
            of the energy band which the flux was computed in.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only
            be specified if you have a reason to generate the same
            set of random numbers, such as for a test. Default is None,
            which sets the seed based on the system time.
        """
        prng = parse_prng(prng)
        energy = np.asarray(events["energy"])
        if energy.size == 0:
            return events
        earea = self.interpolate_area(energy).value
        idxs = np.logical_and(energy >= refband[0], energy <= refband[1])
        rate = flux / (energy[idxs].sum() * erg_per_keV) * earea[idxs].sum()
        n_ph = prng.poisson(lam=rate * exp_time)
        fak = float(n_ph) / energy.size
        if fak > 1.0:
            mylog.error(
                "Number of events in sample: %d, Number of events " "wanted: %d",
                n_ph,
                energy.size,
            )
            raise ValueError(
                "This combination of exposure time and effective "
                "area will result in more photons being drawn "
                "than are available in the sample!!!"
            )
        w = earea / self.max_area
        randvec = prng.uniform(size=energy.size)
        eidxs = prng.permutation(np.where(randvec < w)[0])[:n_ph].astype("int64")
        mylog.debug("%d events detected out of %d.", n_ph, energy.size)
        for key in events:
            events[key] = np.asarray(events[key][eidxs])
        return events

    def plot(
        self,
        xscale="log",
        yscale="log",
        xlabel=None,
        ylabel=None,
        fig=None,
        ax=None,
        **kwargs,
    ):
        """
        Make a quick plot of the effective area curve.

        Parameters
        ----------
        xscale : string
            The scale of the x-axis. "linear" or "log".
        yscale : string
            The scale of the y-axis. "linear" or "log".
        xlabel : string
            The label of the x-axis. Default: "E (keV)"
        ylabel : string
            The label of the y-axis. Default: "$\\mathrm{A\\ (cm^2)}$"
        fig : :class:`~matplotlib.figure.Figure`, optional
            The figure to place the plot in. If not supplied,
            one will be created.
        ax : :class:`~matplotlib.axes.Axes`, optional
            The axes to place the plot in. If not supplied,
            one will be created.

        Notes
        -----
        All other arguments are passed to the call to
        :meth:`~matplotlib.axes.Axes.plot`.

        Returns
        -------

        A tuple of the :class:`~matplotlib.figure.Figure` and
        :class:`~matplotlib.axes.Axes` objects.
        """
        import matplotlib.pyplot as plt

        if xlabel is None:
            xlabel = "E (keV)"
        if ylabel is None:
            ylabel = "$\\mathrm{A\\ (cm^2)}$"
        if fig is None:
            fig = plt.figure(figsize=(10, 10))
        if ax is None:
            ax = fig.add_subplot(111)
        ax.plot(self.emid, self.eff_area, **kwargs)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return fig, ax


class FlatResponse(AuxiliaryResponseFile):
    """
    A flat effective area response.

    Parameters
    ----------
    emin : float
        The minimum energy of the response in keV.
    emax : float
        The maximum energy of the response in keV.
    area : float
        The effective area in cm**2.
    nbins : integer
        The number of bins in the response file.

    Examples
    --------
    >>> arf = FlatResponse(0.1, 10.0, 3000.0, 10000)
    """

    def __init__(self, emin, emax, area, nbins):
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        area = parse_value(area, "cm**2")
        self.filename = "flat_response"
        self.de = (emax - emin) / nbins
        self.ebins = np.linspace(emin, emax, nbins + 1)
        self.elo = self.ebins[:-1]
        self.ehi = self.ebins[1:]
        self.emid = 0.5 * (self.elo + self.ehi)
        self.eff_area = area * np.ones(nbins)
        self.max_area = area

    def interpolate_area(self, energy):
        return u.Quantity(self.max_area * np.ones_like(energy), "cm**2")


class RedistributionMatrixFile:
    r"""
    A class for redistribution matrix files (RMFs).

    Parameters
    ----------
    filename : string
        The filename of the RMF to be read.

    Examples
    --------
    >>> rmf = RedistributionMatrixFile("xrs_hdxi.rmf")
    """

    def __init__(self, filename):
        self.filename = get_data_file(filename)
        self.handle = fits.open(self.filename, memmap=True)
        if "MATRIX" in self.handle:
            self.mat_key = "MATRIX"
        elif "SPECRESP MATRIX" in self.handle:
            self.mat_key = "SPECRESP MATRIX"
        else:
            raise RuntimeError(
                f"Cannot find the response matrix in the RMF "
                f"file {filename}! It should be named "
                f'"MATRIX" or "SPECRESP MATRIX".'
            )
        self.header = self.handle[self.mat_key].header
        self.num_mat_columns = len(self.handle[self.mat_key].columns)
        self.ebounds_header = self.handle["EBOUNDS"].header
        self.weights = np.array([np.nansum(w) for w in self.data["MATRIX"]])
        self.elo = self.data["ENERG_LO"]
        self.ehi = self.data["ENERG_HI"]
        self.ebins = np.append(self.data["ENERG_LO"], self.data["ENERG_HI"][-1])
        self.emid = 0.5 * (self.elo + self.ehi)
        self.de = self.ehi - self.elo
        self.n_e = self.elo.size
        self.n_ch = self.header["DETCHANS"]
        num = 0
        for i in range(1, self.num_mat_columns + 1):
            if self.header[f"TTYPE{i}"] == "F_CHAN":
                num = i
                break
        self.cmin = self.header.get(f"TLMIN{num}", 1)
        self.cmax = self.header.get(f"TLMAX{num}", self.n_ch)

    @classmethod
    def from_instrument(cls, name):
        """
        Return an :class:`~soxs.instrument.RedistributionMatrixFile`
        object from the name of an existing instrument
        specification in SOXS.

        Parameters
        ----------
        name : string
            The name of the instrument specification to use
            to obtain the RMF object from.

        Examples
        --------
        >>> arf = soxs.RedistributionMatrixFile.from_instrument("hdxi")
        """
        instr = instrument_registry.get(name, None)
        if instr is None:
            raise KeyError(f"Instrument '{name}' not in registry!")
        return cls(instr["rmf"])

    @property
    def chan_type(self):
        if "CHANTYPE" in self.header:
            ctype = self.header["CHANTYPE"]
        elif "CHANTYPE" in self.ebounds_header:
            ctype = self.ebounds_header["CHANTYPE"]
        else:
            raise KeyError("'CHANTYPE' not specified in RMF!!")
        return ctype.upper()

    @property
    def data(self):
        return self.handle[self.mat_key].data

    @property
    def ebounds_data(self):
        return self.handle["EBOUNDS"].data

    def __str__(self):
        return self.filename

    def _make_channels(self, k):
        # build channel number list associated to array value,
        # there are groups of channels in rmfs with nonzero probabilities
        true_channel = []
        f_chan = ensure_numpy_array(np.nan_to_num(self.data["F_CHAN"][k]))
        n_chan = ensure_numpy_array(np.nan_to_num(self.data["N_CHAN"][k]))
        for start, nchan in zip(f_chan, n_chan):
            if nchan == 0:
                true_channel.append(start)
            else:
                true_channel += list(range(start, start + nchan))
        return np.array(true_channel)

    def eb_to_ch(self, energy):
        energy = parse_value(energy, "keV")
        idxs = np.searchsorted(self.ebounds_data["E_MIN"], energy) - 1
        return np.arange(self.n_ch)[idxs] + self.cmin

    def ch_to_eb(self, channels, prng=None):
        prng = parse_prng(prng)
        emin = self.ebounds_data["E_MIN"]
        emax = self.ebounds_data["E_MAX"]
        de = emax - emin
        ch = channels - self.cmin
        e = emin[ch] + prng.uniform(size=channels.size) * de[ch]
        return e

    def scatter_energies(self, events, prng=None):
        """
        Scatter photon energies with the RMF and produce the
        corresponding channel values.

        Parameters
        ----------
        events : dict of np.ndarrays
            The energies and positions of the photons.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only
            be specified if you have a reason to generate the same
            set of random numbers, such as for a test. Default is None,
            which sets the seed based on the system time.
        """
        prng = parse_prng(prng)
        eidxs = np.argsort(events["energy"])
        sorted_e = np.asarray(events["energy"])[eidxs]

        detected_channels = []

        # run through all photon energies and find which bin they go in
        fcurr = 0
        last = sorted_e.shape[0]

        emin = sorted_e[0]
        emax = sorted_e[-1]

        pbar = tqdm(leave=True, total=last, desc="Scattering energies ")
        for (k, low), high in zip(enumerate(self.elo), self.ehi):
            if high < emin or low > emax:
                continue
            e = sorted_e[fcurr:last]
            nn = np.logical_and(low <= e, e < high).sum()
            if nn == 0:
                continue
            # weight function for probabilities from RMF
            weights = np.nan_to_num(np.float64(self.data["MATRIX"][k]))
            weights /= weights.sum()
            true_channel = self._make_channels(k)
            if len(true_channel) > 0:
                channel_ind = prng.choice(len(weights), size=nn, p=weights)
                detected_channels.append(true_channel[channel_ind])
                fcurr += nn
                pbar.update(nn)

        pbar.close()

        for key in events:
            events[key] = events[key][eidxs]
        events[self.chan_type] = np.concatenate(detected_channels)
        events["soxs_energy"] = events["energy"].copy()
        events["energy"] = self.ch_to_eb(events[self.chan_type], prng=prng)

        return events

    def convolve_spectrum(self, cspec, exp_time, noisy=True, prng=None, rate=False):
        from soxs.spectra import ConvolvedSpectrum

        prng = parse_prng(prng)
        exp_time = parse_value(exp_time, "s")
        if isinstance(cspec, ConvolvedSpectrum):
            counts = cspec.flux.value * exp_time * cspec.de.value
            if (
                len(cspec.emid) == self.n_e
                and np.isclose(cspec.ebins.value, self.ebins).all()
            ):
                spec = counts
            else:
                spec = regrid_spectrum(self.ebins, cspec.ebins.value, counts)
        else:
            spec = np.asarray(cspec) * exp_time
        conv_spec = np.zeros(self.n_ch)
        pbar = tqdm(leave=True, total=self.n_e, desc="Convolving spectrum ")
        if not isinstance(self.data["MATRIX"], fits.column._VLF) and np.all(
            self.data["N_GRP"] == 1
        ):
            # We can do things a bit faster if there is only one group each
            f_chan = ensure_numpy_array(np.nan_to_num(self.data["F_CHAN"])) - self.cmin
            n_chan = ensure_numpy_array(np.nan_to_num(self.data["N_CHAN"]))
            mat = np.nan_to_num(np.float64(self.data["MATRIX"]))
            for k in range(self.n_e):
                weights = mat[k, :]
                weights /= weights.sum()
                conv_spec[f_chan[k] : f_chan[k] + n_chan[k]] += (
                    spec[k] * weights[: n_chan[k]]
                )
                pbar.update()
        else:
            # Otherwise, we have to go step-by-step
            for k in range(self.n_e):
                f_chan = (
                    ensure_numpy_array(np.nan_to_num(self.data["F_CHAN"][k]))
                    - self.cmin
                )
                n_chan = ensure_numpy_array(np.nan_to_num(self.data["N_CHAN"][k]))
                weights = np.nan_to_num(np.float64(self.data["MATRIX"][k]))
                weights /= weights.sum()
                f1 = 0
                for n, f in zip(n_chan, f_chan):
                    conv_spec[f : f + n] += spec[k] * weights[f1 : f1 + n]
                    f1 += n
                pbar.update()
        pbar.close()
        if noisy:
            conv_spec = prng.poisson(lam=conv_spec)
        if rate:
            conv_spec /= exp_time
        return conv_spec
