import copy
import os
import shutil
import subprocess
import tempfile
import warnings

import astropy.units as u
import h5py
import numpy as np
from astropy.modeling.functional_models import Gaussian1D
from astropy.table import QTable
from pathlib import Path, PurePath

from soxs.constants import ckms, erg_per_keV, hc, sigma_to_fwhm, sqrt2pi
from soxs.spectra.foreground_absorption import tbabs_cross_section, wabs_cross_section
from soxs.utils import (
    issue_deprecation_warning,
    line_width_equiv,
    mylog,
    parse_prng,
    parse_value,
    regrid_spectrum,
)


class Energies(u.Quantity):
    def __new__(cls, energy, flux):
        ret = u.Quantity.__new__(cls, energy, unit="keV")
        ret.flux = u.Quantity(flux, "erg/(cm**2*s)")
        return ret


def _generate_energies(spec, t_exp, rate, prng, binscale, quiet=False):
    cumspec = spec.cumspec
    n_ph = prng.poisson(t_exp * rate)
    if not quiet:
        mylog.info("Creating %d energies from this spectrum.", n_ph)
    randvec = prng.uniform(size=n_ph)
    randvec.sort()
    if binscale in ["linear", "custom"]:
        bin_edges = spec.ebins.value
    elif binscale == "log":
        bin_edges = np.log10(spec.ebins.value)
    e = np.interp(randvec, cumspec, bin_edges)
    if binscale == "log":
        e = 10**e
    if not quiet:
        mylog.info("Finished creating energies.")
    return e


def distance_handler(redshift, dist, cosmology):
    from astropy.cosmology import Planck18

    if cosmology is None:
        cosmology = Planck18
    if redshift == 0.0 and dist is None:
        raise ValueError("Either 'redshift' must be > 0 or 'dist' cannot be None!")
    if dist is None:
        D_A = cosmology.angular_diameter_distance(redshift).to_value("cm")
    else:
        D_A = u.Quantity(parse_value(dist, "kpc"), "kpc").to_value("cm")
    return D_A


class BaseSpectrum:
    _units = ""  # set by subclasses

    def __init__(self, ebins, spec, binscale="linear"):
        self.ebins = u.Quantity(ebins, "keV")
        self.emid = 0.5 * (self.ebins[1:] + self.ebins[:-1])
        self._spec = u.Quantity(spec, self._units)
        self.nbins = len(self.emid)
        self.de = np.diff(self.ebins)
        self.binscale = binscale
        self._compute_totals()
        self.noisy = False
        self._spec_err = None
        self.exp_time = None
        self.arf = None
        self.rmf = None

    def __copy__(self):
        # Create a new instance and copy attributes
        return type(self)(self.ebins, self._spec.value, binscale=self.binscale)

    def __deepcopy__(self, memo):
        # Create a new instance and deep copy attributes
        return type(self)(
            copy.deepcopy(self.ebins, memo),
            copy.deepcopy(self._spec, memo),
            binscale=self.binscale,
        )

    def _compute_totals(self):
        self._energy_spec = self._spec * self.emid.to("erg") / (1.0 * u.photon)
        self._total_spec = (self._spec * self.de).sum()
        self._total_energy_spec = (self._energy_spec * self.de).sum()
        cumspec = np.cumsum((self._spec * self.de).value)
        cumspec = np.insert(cumspec, 0, 0.0)
        cumspec /= cumspec[-1]
        self.cumspec = cumspec
        self.func = lambda e: np.interp(e, self.emid.value, self._spec.value)
        self._binned_spec = self._spec * self.de
        self._binned_energy_spec = self._energy_spec * self.de

    def _check_binning_units(self, other):
        if self.nbins != other.nbins or not np.isclose(self.ebins.value, other.ebins.value).all():
            raise RuntimeError("Energy binning for these two spectra is not the same!!")
        if self._units != other._units:
            raise RuntimeError("The units for these two spectra are not the same!")

    def _compute_waves(self):
        self._wvbins = self.ebins.to("angstrom", equivalencies=u.spectral())
        self._wvmid = 0.5 * (self._wvbins[1:] + self._wvbins[:-1])
        self._dwv = -np.diff(self._wvbins)

    def _compute_freqs(self):
        self._fbins = self.ebins.to("Hz", equivalencies=u.spectral())
        self._fmid = 0.5 * (self._fbins[1:] + self._fbins[:-1])
        self._df = np.diff(self._fbins)

    _wvbins = None

    @property
    def wvbins(self):
        if self._wvbins is None:
            self._compute_waves()
        return self._wvbins

    _wvmid = None

    @property
    def wvmid(self):
        if self._wvmid is None:
            self._compute_waves()
        return self._wvmid

    _dwv = None

    @property
    def dwv(self):
        if self._dwv is None:
            self._compute_waves()
        return self._dwv

    _fbins = None

    @property
    def fbins(self):
        if self._fbins is None:
            self._compute_freqs()
        return self._fbins

    _fmid = None

    @property
    def fmid(self):
        if self._fmid is None:
            self._compute_freqs()
        return self._fmid

    _df = None

    @property
    def df(self):
        if self._df is None:
            self._compute_freqs()
        return self._df

    def __add__(self, other):
        self._check_binning_units(other)
        return type(self)(self.ebins, self._spec + other._spec, binscale=self.binscale)

    def __sub__(self, other):
        self._check_binning_units(other)
        return type(self)(self.ebins, self._spec - other._spec, binscale=self.binscale)

    def __mul__(self, other):
        if hasattr(other, "eff_area"):
            return ConvolvedSpectrum.convolve(self, other)
        else:
            return type(self)(self.ebins, other * self._spec, binscale=self.binscale)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return type(self)(self.ebins, self._spec / other)

    __div__ = __truediv__

    def __repr__(self):
        s = f"{type(self).__name__} ({self.ebins[0]} - {self.ebins[-1]})\n"
        s += f"    Total Values:\n    {self._total_spec}\n    {self._total_energy_spec}\n"
        return s

    def __call__(self, e):
        if hasattr(e, "to_astropy"):
            e = e.to_astropy()
        if isinstance(e, u.Quantity):
            e = e.to("keV").value
        return u.Quantity(self.func(e), self._units)

    def regrid_spectrum(self, emin, emax, nbins, binscale="linear", vlos=0.0, vtot=None):
        """
        Regrid an existing spectrum to a new energy binning and
        return a new spectrum.

        Parameters
        ----------
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy in the band, in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy in the band, in keV.
        nbins : integer
            The number of bins in the spectrum.
        binscale : string, optional
            The scale of the energy binning: "linear" or "log".
            Default: "linear"
        vlos : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The line-of-sight velocity of the source, in km/s.
            Used for Doppler-shifting the spectrum. Default: 0.0
        vtot : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The total velocity magnitude of the source, in km/s.
            If not set, it is assumed to be equal to the absolute
            value of the vlos argument. This value is mainly relevant
            if the velocity is relativistic and is not completely
            along the sight line. Default: None
        """
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        vlos = parse_value(vlos, "km/s") / ckms
        if vtot is None:
            vtot = np.abs(vlos)
        else:
            vtot = parse_value(vtot, "km/s") / ckms
        shift = np.sqrt(1.0 - vtot**2) / (1.0 + vlos)
        if binscale == "linear":
            ebins = np.linspace(emin, emax, nbins + 1)
        elif binscale == "log":
            ebins = np.logspace(np.log10(emin), np.log10(emax), nbins + 1)
        de = np.diff(ebins)
        spec_new = (
            shift
            * shift
            * shift
            * regrid_spectrum(ebins / shift, self.ebins.value, self._spec.value * self.de.value)
            / de
        )
        return type(self)(ebins, spec_new, binscale=binscale)

    def restrict_within_band(self, emin=None, emax=None):
        """
        Zeros out the spectrum outside a specified energy band.

        Parameters
        ----------
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The minimum energy in the band, in keV. Default is to use the lowest energy.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The maximum energy in the band, in keV. Default is to use the highest energy.
        """
        if emin is not None:
            emin = parse_value(emin, "keV")
            self._spec[self.emid.value < emin] = 0.0
        if emax is not None:
            emax = parse_value(emax, "keV")
            self._spec[self.emid.value > emax] = 0.0
        self._compute_totals()

    def _get_total_in_band(self, emin, emax):
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        range = np.logical_and(self.emid.value >= emin, self.emid.value <= emax)
        pflux = self._binned_spec[range].sum()
        eflux = self._binned_energy_spec[range].sum()
        return pflux, eflux

    @classmethod
    def _from_ext_model(cls, ebins, flux):
        ebins = np.array(ebins)
        de = np.diff(ebins)
        dloge = np.diff(np.log10(ebins))
        flux = np.array(flux) / de
        if np.isclose(de, de[0]).all():
            binscale = "linear"
        elif np.isclose(dloge, dloge[0]).all():
            binscale = "log"
        else:
            binscale = "custom"
        return cls(ebins, flux, binscale=binscale)

    @classmethod
    def from_constant(cls, const_value, emin, emax, nbins, binscale="linear"):
        """
        Create a spectrum from a constant model.

        Parameters
        ----------
        const_value : float
            The value of the constant in the units of the spectrum.
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy of the spectrum in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy of the spectrum in keV.
        nbins : integer
            The number of bins in the spectrum.
        binscale : string, optional
            The scale of the energy binning: "linear" or "log".
            Default: "linear"
        """
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        if binscale == "linear":
            ebins = np.linspace(emin, emax, nbins + 1)
        elif binscale == "log":
            ebins = np.logspace(np.log10(emin), np.log10(emax), nbins + 1)
        flux = const_value * np.ones(nbins)
        return cls(ebins, flux, binscale=binscale)

    @classmethod
    def from_powerlaw(cls, photon_index, redshift, norm, emin, emax, nbins, binscale="linear"):
        """
        Create a spectrum from a power-law model. The form of the model
        is F(E) = norm*(e*(1+z)/(1 keV))**-photon_index.

        Parameters
        ----------
        photon_index : float
            The photon index of the source.
        redshift : float
            The redshift of the source.
        norm : float
            The normalization of the source in units of
            the spectrum, photons/s/cm**2/keV, at 1 keV
            in the source frame.
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy of the spectrum in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy of the spectrum in keV.
        nbins : integer
            The number of bins in the spectrum.
        binscale : string, optional
            The scale of the energy binning: "linear" or "log".
            Default: "linear"
        """
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        if binscale == "linear":
            ebins = np.linspace(emin, emax, nbins + 1)
        elif binscale == "log":
            ebins = np.logspace(np.log10(emin), np.log10(emax), nbins + 1)
        emid = 0.5 * (ebins[1:] + ebins[:-1])
        flux = norm * (emid * (1.0 + redshift)) ** (-photon_index)
        return cls(ebins, flux, binscale=binscale)

    @classmethod
    def from_file(cls, filename):
        """
        Read a spectrum from an ASCII ECSV, FITS, or HDF5 file, in
        the forms written by ``write_ascii_file``, ``write_fits_file``,
        and ``write_h5_file``.

        Parameters
        ----------
        filename : string
            The path to the file containing the spectrum.
        """
        kwargs = {}
        try:
            # First try reading the file as HDF5
            with h5py.File(filename, "r") as f:
                flux = f["spectrum"][()]
                nbins = flux.size
                binscale = f.attrs.get("binscale", None)
                if binscale is None:
                    warnings.warn(
                        "This spectrum file was written without a 'binscale' attribute! "
                        "A linear energy binning will be assumed.",
                        stacklevel=2,
                    )
                    binscale = "linear"
                if binscale == "linear":
                    ebins = np.linspace(f["emin"][()], f["emax"][()], nbins + 1)
                elif binscale == "log":
                    ebins = np.logspace(np.log10(f["emin"][()]), np.log10(f["emax"][()]), nbins + 1)
                for key in ["arf", "rmf", "noisy", "exp_time"]:
                    kwargs[key] = f.attrs.get(key, None)
        except OSError:
            # If that fails, try reading the file as ASCII or FITS
            # using AstroPy QTable
            try:
                t = QTable.read(filename, format="fits")
            except OSError:
                t = QTable.read(filename, format="ascii.ecsv")
            ebins = np.append(t["elo"].value, t["ehi"].value[-1])
            flux = t["flux"].value
            binscale = t.meta.get("binscale", t.meta.get("BINSCALE", None))
            if binscale is None:
                warnings.warn(
                    "This spectrum file was written without a 'binscale' attribute! "
                    "A linear energy binning will be assumed.",
                    stacklevel=2,
                )
                binscale = "linear"
            units = t.meta.get("units", t.meta.get("UNITS", None))
            if units is None or units != cls._units:
                raise ValueError(
                    f"Spectrum units in file ({t.meta['UNITS']}) do not match "
                    f"the expected units ({cls._units})!"
                ) from None
            for key in ["arf", "rmf", "noisy", "exp_time"]:
                kwargs[key] = t.meta.get(key, t.meta.get(key.upper(), None))
        if kwargs["arf"] is not None:
            arf = kwargs.pop("arf")
            return cls(
                ebins,
                flux,
                arf,
                binscale=binscale,
                **kwargs,
            )
        else:
            return cls(ebins, flux, binscale=binscale)

    def _new_spec_from_band(self, emin, emax):
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        band = np.logical_and(self.ebins.value >= emin, self.ebins.value <= emax)
        idxs = np.where(band)[0]
        ebins = self.ebins.value[idxs]
        flux = self._spec.value[idxs[:-1]]
        return ebins, flux

    def new_spec_from_band(self, emin, emax):
        """
        Create a new :class:`~soxs.spectra.base.Spectrum` object
        from a subset of an existing one defined by a particular
        energy band.

        Parameters
        ----------
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy of the band in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy of the band in keV.
        """
        ebins, flux = self._new_spec_from_band(emin, emax)
        return type(self)(ebins, flux, binscale=self.binscale)

    def _write_fits_or_ascii(self):
        meta = {"binscale": self.binscale, "noisy": self.noisy, "units": self._units}
        if self.arf is not None:
            meta["arf"] = self.arf.filename
        if self.rmf is not None:
            meta["rmf"] = self.rmf.filename
        if self.exp_time is not None:
            meta["exp_time"] = self.exp_time.value
        t = QTable(
            [self.ebins[:-1], self.ebins[1:], self._spec],
            names=["elo", "ehi", "flux"],
            meta=meta,
        )
        return t

    def write_ascii_file(self, specfile, overwrite=False):
        """
        Write the spectrum to an ASCII file in the ECSV format
        (https://docs.astropy.org/en/stable/io/ascii/ecsv.html).

        Parameters
        ----------
        specfile : string
            The filename to write the file to.
        overwrite : boolean, optional
            Whether to overwrite an existing
            file with the same name. Default: False
        """
        t = self._write_fits_or_ascii()
        t.write(specfile, overwrite=overwrite, format="ascii.ecsv")

    def write_file(self, specfile, overwrite=False):
        """
        Write a :class:`~soxs.Spectrum` object to disk, in
        any of three formats, which will be determined by
        the chosen suffix:
        ASCII ECSV: .dat, .txt. or .ecsv
        HDF5: .hdf5 or .h5
        FITS: .fits

        Parameters
        ----------
        specfile : string
            The name of the file to write to.
        overwrite : boolean, optional
            Whether to overwrite an existing
            file with the same name. Default: False
        """
        suffix = PurePath(specfile).suffix.lower()
        if suffix in [".dat", ".txt", ".ecsv"]:
            mylog.info("Writing ASCII ECSV file %s.", specfile)
            self.write_ascii_file(specfile, overwrite=overwrite)
        elif suffix in [".hdf5", ".h5"]:
            mylog.info("Writing HDF5 file %s.", specfile)
            self.write_hdf5_file(specfile, overwrite=overwrite)
        elif suffix == ".fits":
            mylog.info("Writing FITS file %s.", specfile)
            self.write_fits_file(specfile, overwrite=overwrite)
        else:
            raise NotImplementedError(f"Unknown file suffix {suffix}!")

    def write_hdf5_file(self, specfile, overwrite=False):
        """
        Write the spectrum to an HDF5 file.

        Parameters
        ----------
        specfile : string
            The filename to write the file to.
        overwrite : boolean, optional
            Whether to overwrite an existing
            file with the same name. Default: False
        """
        if Path(specfile).exists() and not overwrite:
            raise OSError(f"File {specfile} exists and overwrite=False!")
        with h5py.File(specfile, "w") as f:
            f.create_dataset("emin", data=self.ebins[0].value)
            f.create_dataset("emax", data=self.ebins[-1].value)
            f.create_dataset("spectrum", data=self._spec.value)
            f.attrs["binscale"] = self.binscale
            f.attrs["noisy"] = self.noisy
            f.attrs["units"] = self._units
            if self.arf is not None:
                f.attrs["arf"] = self.arf.filename
            if self.rmf is not None:
                f.attrs["rmf"] = self.rmf.filename
            if self.exp_time is not None:
                f.attrs["exp_time"] = self.exp_time.value

    def write_fits_file(self, specfile, overwrite=False):
        """
        Write the spectrum to a FITS file.

        Parameters
        ----------
        specfile : string
            The filename to write the file to.
        overwrite : boolean, optional
            Whether to overwrite an existing
            file with the same name. Default: False
        """
        t = self._write_fits_or_ascii()
        t.write(specfile, overwrite=overwrite, format="fits")

    def add_emission_line(self, line_center, line_width, line_amp, line_type="gaussian"):
        """
        Add an emission line to this spectrum.

        Parameters
        ----------
        line_center : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The line center position in units of keV, in the observer frame.
        line_width : one or more float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The line width (FWHM) in units of keV, in the observer frame. Can also
            input the line width in units of velocity in the rest frame. For the Voigt
            profile, a list, tuple, or array of two values should be provided since there
            are two line widths, the Lorentzian and the Gaussian (in that order).
        line_amp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The integrated line amplitude in the units of the spectrum
        line_type : string, optional
            The line profile type. Default: "gaussian"
        """
        line_center = parse_value(line_center, "keV")
        line_width = parse_value(line_width, "keV", equivalence=line_width_equiv(line_center))
        line_amp = parse_value(line_amp, self._units)
        if line_type == "gaussian":
            sigma = line_width / sigma_to_fwhm
            line_amp /= sqrt2pi * sigma
            f = Gaussian1D(line_amp, line_center, sigma)
        else:
            raise NotImplementedError(f'Line profile type "{line_type}" not implemented!')
        self._spec += u.Quantity(f(self.emid.value), self._units)
        self._compute_totals()

    def add_absorption_line(self, line_center, line_width, equiv_width, line_type="gaussian"):
        """
        Add an absorption line to this spectrum.

        Parameters
        ----------
        line_center : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The line center position in units of keV, in the observer frame.
        line_width : one or more float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The line width (FWHM) in units of keV, in the observer frame. Can also
            input the line width in units of velocity in the rest frame. For the Voigt
            profile, a list, tuple, or array of two values should be provided since there
            are two line widths, the Lorentzian and the Gaussian (in that order).
        equiv_width : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The equivalent width of the line, in units of milli-Angstrom
        line_type : string, optional
            The line profile type. Default: "gaussian"
        """
        line_center = parse_value(line_center, "keV")
        line_width = parse_value(line_width, "keV", equivalence=line_width_equiv(line_center))
        equiv_width = parse_value(equiv_width, "1.0e-3*angstrom")  # in milliangstroms
        equiv_width *= 1.0e-3  # convert to angstroms
        if line_type == "gaussian":
            sigma = line_width / sigma_to_fwhm
            B = equiv_width * line_center * line_center
            B /= hc * sqrt2pi * sigma
            f = Gaussian1D(B, line_center, sigma)
        else:
            raise NotImplementedError(f"Line profile type '{line_type}' not implemented!")
        self._spec *= np.exp(-f(self.emid.value))
        self._compute_totals()

    def apply_foreground_absorption(self, nH, model="wabs", redshift=0.0, abund_table="angr"):
        """
        Given a hydrogen column density, apply
        galactic foreground absorption to the spectrum.

        Parameters
        ----------
        nH : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The hydrogen column in units of 10**22 atoms/cm**2
        model : string, optional
            The model for absorption to use. Options are "wabs"
            (Wisconsin, Morrison and McCammon; ApJ 270, 119) or
            "tbabs" (Tuebingen-Boulder, Wilms, J., Allen, A., &
            McCray, R. 2000, ApJ, 542, 914). Default: "wabs".
        redshift : float, optional
            The redshift of the absorbing material. Default: 0.0
        abund_table : string
            The abundance table to be used for abundances in the
            absorption model. Default is set in the SOXS
            configuration file, the default for which is "angr".
            Built-in options are:
            "angr" : from Anders E. & Grevesse N. (1989, Geochimica et
            Cosmochimica Acta 53, 197)
            "aspl" : from Asplund M., Grevesse N., Sauval A.J. & Scott
            P. (2009, ARAA, 47, 481)
            "feld" : from Feldman U. (1992, Physica Scripta, 46, 202)
            "wilm" : from Wilms, Allen & McCray (2000, ApJ 542, 914
            except for elements not listed which are given zero abundance)
            "lodd" : from Lodders, K (2003, ApJ 591, 1220)
            "cl17.03" : the default abundance table in Cloudy 17.03
        """
        nH = parse_value(nH, "1.0e22*cm**-2")
        e = self.emid.value * (1.0 + redshift)
        if model == "wabs":
            sigma = wabs_cross_section(e)
        elif model == "tbabs":
            if not isinstance(abund_table, str):
                raise ValueError(
                    "Must supply a string corresponding to one of "
                    "the abundance tables provided in SOXS for "
                    "the TBabs model!"
                )
            sigma = tbabs_cross_section(e, abund_table=abund_table)
        self._spec *= np.exp(-nH * 1.0e22 * sigma)
        self._compute_totals()

    def plot(
        self,
        lw=2,
        xmin=None,
        xmax=None,
        ymin=None,
        ymax=None,
        xscale=None,
        yscale=None,
        label=None,
        fontsize=18,
        fig=None,
        ax=None,
        **kwargs,
    ):
        """
        Make a quick Matplotlib plot of the spectrum. A Matplotlib
        figure and axis is returned.

        Parameters
        ----------
        lw : float, optional
            The width of the lines in the plots. Default: 2.0 px.
        xmin : float, optional
            The left-most energy in keV to plot. Default is the
            minimum value in the spectrum.
        xmax : float, optional
            The right-most energy in keV to plot. Default is the
            maximum value in the spectrum.
        ymin : float, optional
            The lower extent of the y-axis. By default, it is set automatically.
        ymax : float, optional
            The upper extent of the y-axis. By default, it is set automatically.
        xscale : string, optional
            The scaling of the x-axis of the plot. Default: "log"
        yscale : string, optional
            The scaling of the y-axis of the plot. Default: "log"
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

        Returns
        -------
        A tuple of the :class:`~matplotlib.figure.Figure` and the
        :class:`~matplotlib.axes.Axes` objects.
        """
        import matplotlib.pyplot as plt

        if fig is None:
            fig = plt.figure(figsize=(10, 10))
        if xscale is None:
            if ax is None:
                xscale = "log"
            else:
                xscale = ax.get_xscale()
        if yscale is None:
            if ax is None:
                yscale = "log"
            else:
                yscale = ax.get_yscale()
        if ax is None:
            ax = fig.add_subplot(111)
        if self.noisy:
            ax.errorbar(
                self.emid,
                self._spec,
                xerr=0.5 * self.de,
                yerr=self._spec_err,
                fmt=".",
                lw=lw,
                ms=lw,
                label=label,
                **kwargs,
            )
        else:
            ax.plot(self.emid, self._spec, lw=lw, label=label, **kwargs)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel("Energy (keV)", fontsize=fontsize)
        yunit = u.Unit(self._units).to_string("latex").replace("{}^{\\prime}", "arcmin")
        ax.set_ylabel(f"Spectrum ({yunit})", fontsize=fontsize)
        ax.tick_params(reset=True, axis="both", labelsize=fontsize)
        return fig, ax


class Spectrum(BaseSpectrum):
    _units = "photon/(cm**2*s*keV)"

    @property
    def flux(self):
        return self._spec

    @property
    def energy_flux(self):
        return self._energy_spec

    @property
    def total_flux(self):
        return self._total_spec

    @property
    def total_energy_flux(self):
        return self._total_energy_spec

    @property
    def binned_flux(self):
        return self._binned_spec

    @property
    def binned_energy_flux(self):
        return self._binned_energy_spec

    @property
    def flux_per_wavelength(self):
        return self.binned_flux / self.dwv

    @property
    def energy_flux_per_wavelength(self):
        return self.binned_energy_flux / self.dwv

    @property
    def flux_per_frequency(self):
        return self.binned_flux / self.df

    @property
    def energy_flux_per_frequency(self):
        return self.binned_energy_flux / self.df

    def rescale_flux(self, new_flux, emin=None, emax=None, flux_type="photons"):
        """
        Rescale the flux of the spectrum, optionally using
        a specific energy band.

        Parameters
        ----------
        new_flux : float
            The new flux in units of photons/s/cm**2.
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The minimum energy of the band to consider,
            in keV. Default: Use the minimum energy of
            the entire spectrum.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The maximum energy of the band to consider,
            in keV. Default: Use the maximum energy of
            the entire spectrum.
        flux_type : string, optional
            The units of the flux to use in the rescaling:
                "photons": photons/s/cm**2
                "energy": erg/s/cm**2
        """
        if emin is None:
            emin = self.ebins[0].value
        if emax is None:
            emax = self.ebins[-1].value
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        idxs = np.logical_and(self.emid.value >= emin, self.emid.value <= emax)
        if flux_type == "photons":
            f = self._binned_spec[idxs].sum()
        elif flux_type == "energy":
            f = self._binned_energy_spec[idxs].sum()
        self._spec *= new_flux / f.value
        self._compute_totals()

    def get_flux_in_band(self, emin, emax):
        """
        Determine the total flux within a band specified
        by an energy range.

        Parameters
        ----------
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy in the band, in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy in the band, in keV.

        Returns
        -------
        A tuple of values for the flux/intensity in the
        band: the first value is in terms of the photon
        rate, the second value is in terms of the energy rate.
        """
        return self._get_total_in_band(emin, emax)

    def get_lum_in_band(self, emin, emax, redshift=0.0, dist=None, cosmology=None):
        """
        Determine the total photon count rate and luminosity
        within a band specified by an energy range. Either a redshift or
        or an explicit distance must be specified in order to convert
        from the observer frame to the source frame.

        Parameters
        ----------
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy in the band, in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy in the band, in keV.
        redshift : float
            The redshift to the source, assuming it is cosmological.
        dist : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The distance to the source, if it is not cosmological. If a unit
            is not specified, it is assumed to be in kpc.
        cosmology : :class:`~astropy.cosmology.Cosmology` object
            An AstroPy cosmology object used to determine the luminosity
            distance if needed. If not set, the default is the Planck 2018
            cosmology.

        Returns
        -------
        A tuple of values for the luminosity in the band: the first
        value is in terms of the photon rate, the second value is
        in terms of the energy rate.
        """
        D_L = distance_handler(redshift, dist, cosmology) * (1.0 + redshift) ** 2
        pflux, eflux = self.get_flux_in_band(emin / (1.0 + redshift), emax / (1.0 + redshift))
        dist_fac = 4.0 * np.pi * D_L**2
        elum = dist_fac * eflux
        plum = dist_fac * pflux / (1.0 + redshift)
        return plum, elum

    def generate_energies(self, t_exp, area, prng=None, quiet=False):
        """
        Generate photon energies from this spectrum
        given an exposure time and effective area.

        Parameters
        ----------
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The effective area in cm**2. If one is creating
            events for a SIMPUT file, a constant should be
            used, and it must be large enough so that a
            sufficiently large sample is drawn for the ARF.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. This will typically only
            be specified if you have a reason to generate the same
            set of random numbers, such as for a test. Default is None,
            which sets the seed based on the system time.
        quiet : boolean, optional
            If True, log messages will not be displayed when
            creating energies. Useful if you have to loop over
            a lot of spectra. Default: False
        """
        t_exp = parse_value(t_exp, "s")
        area = parse_value(area, "cm**2")
        prng = parse_prng(prng)
        rate = area * self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng, self.binscale, quiet=quiet)
        flux = np.sum(energy) * erg_per_keV / t_exp / area
        energies = Energies(energy, flux)
        return energies

    @classmethod
    def from_count_rate_spectrum(cls, spec, redshift=0.0, dist=None, cosmology=None):
        """
        Create a Spectrum from a CountRateSpectrum object, assuming that the latter
        is in the rest frame of the source. To get the Spectrum in the observer
        frame, you must either provide a redshift or a distance to the source to
        convert to flux. Note that the energy bins of the original CountRateSpectrum
        object will be cosmologically redshifted if the redshift is nonzero.

        Parameters
        ----------
        spec : :class:`~soxs.spectra.base.CountRateSpectrum`
            The CountRateSpectrum object to create the Spectrum from.
        redshift : float
            The redshift to the source, assuming it is cosmological.
        dist : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The distance to the source, if it is not cosmological. If a unit
            is not specified, it is assumed to be in kpc.
        cosmology : :class:`~astropy.cosmology.Cosmology` object
            An AstroPy cosmology object used to determine the luminosity
            distance if needed. If not set, the default is the Planck 2018
            cosmology.
        """
        D_A = distance_handler(redshift, dist, cosmology)
        flux = spec.rate.value / (4.0 * np.pi * D_A**2 * (1.0 + redshift) ** 2)
        return cls(spec.ebins.value / (1.0 + redshift), flux, binscale=spec.binscale)

    @classmethod
    def from_spex_model(cls, spex_session, isect=1):
        """
        Create a spectrum from a SPEX session.

        Parameters
        ----------
        spex_session : :class:`~pyspex.spex.Session`
            The SPEX session object.
        isect : integer, optional
            The sector in the SPEX session to use. Default: 1.
        """
        ss = spex_session.mod_spectrum
        ss.get(isect)
        ehi = ss.energy_upper.to_value("keV")
        elo = ehi - ss.energy_width.to_value("keV")
        ebins = np.append(elo, ehi[-1])
        flux = ss.spectrum.to_value("ph/(bin*s*cm**2)")
        return cls._from_ext_model(ebins, flux)

    @classmethod
    def _from_xspec(cls, xspec_in, emin, emax, nbins, binscale="linear", xspec_settings=None):
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        tmpdir = Path(tempfile.mkdtemp())
        curdir = Path.cwd()
        outscript = []
        if xspec_settings is not None:
            for key, value in xspec_settings.items():
                outscript.append(f"xset {key} {value}\n")
        outscript.extend(xspec_in)
        outscript.append(f"dummyrsp {emin} {emax} {nbins} {binscale[:3]}\n")
        outscript += [
            f"set fp [open {str(tmpdir)}/spec_therm.xspec w+]\n",
            "tclout energies\n",
            "puts $fp $xspec_tclout\n",
            "tclout modval\n",
            "puts $fp $xspec_tclout\n",
            "close $fp\n",
            "quit\n",
        ]
        with open(tmpdir / "xspec.in", "w") as f_xin:
            f_xin.writelines(outscript)
        logfile = curdir / "xspec.log"
        if os.environ["SHELL"] in ["tcsh", "csh"]:
            init_suffix = "csh"
        else:
            init_suffix = "sh"
        sh_file = [
            f"#!{os.environ['SHELL']}\n",
            f". {os.environ['HEADAS']}/headas-init.{init_suffix}\n",
            f"xspec - {str(tmpdir)}/xspec.in\n",
        ]
        with open(tmpdir / "xspec.sh", "w") as f_xs:
            f_xs.writelines(sh_file)
        with open(logfile, "ab") as xsout:
            subprocess.call(
                [os.environ["SHELL"], f"{str(tmpdir)}/xspec.sh"],
                stdout=xsout,
                stderr=xsout,
            )
        with open(tmpdir / "spec_therm.xspec") as f_s:
            lines = f_s.readlines()
        ebins = np.array(lines[0].split()).astype("float64")
        de = np.diff(ebins)
        flux = np.array(lines[1].split()).astype("float64") / de
        shutil.rmtree(str(tmpdir))
        return cls(ebins, flux, binscale=binscale)

    @classmethod
    def from_xspec_script(cls, infile, emin, emax, nbins, binscale="linear", xspec_settings=None):
        """
        Create a model spectrum using a script file as
        input to XSPEC.

        Parameters
        ----------
        infile : string
            Path to the script file to use.
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy of the spectrum in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy of the spectrum in keV.
        nbins : integer
            The number of bins in the spectrum.
        binscale : string, optional
            The scale of the energy binning: "linear" or "log".
            Default: "linear"
        xspec_settings : dict, optional
            A dictionary of XSPEC settings to set before running the
            model, each of which will be run as (e.g.) "xset {key} {value}".
            Default: None.
        """
        with open(infile) as f:
            xspec_in = f.readlines()
        return cls._from_xspec(
            xspec_in,
            emin,
            emax,
            nbins,
            binscale=binscale,
            xspec_settings=xspec_settings,
        )

    @classmethod
    def from_xspec_model(
        cls,
        model_string,
        params,
        emin,
        emax,
        nbins,
        binscale="linear",
        xspec_settings=None,
    ):
        """
        Create a model spectrum using a model string and parameters
        as input to XSPEC.

        Parameters
        ----------
        model_string : string
            The model to create the spectrum from. Use standard XSPEC
            model syntax. Example: "wabs*mekal"
        params : list
            The list of parameters for the model. Must be in the order
            that XSPEC expects.
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy of the spectrum in keV
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy of the spectrum in keV
        nbins : integer
            The number of bins in the spectrum.
        binscale : string, optional
            The scale of the energy binning: "linear" or "log".
            Default: "linear"
        xspec_settings : dict, optional
            A dictionary of XSPEC settings to set before running the
            model, each of which will be run as (e.g.) "xset {key} {value}".
            Default: None.
        """
        xspec_in = []
        model_str = f"{model_string} &"
        for param in params:
            model_str += f" {param} &"
        model_str += " /*"
        xspec_in.append(f"model {model_str}\n")
        return cls._from_xspec(
            xspec_in,
            emin,
            emax,
            nbins,
            binscale=binscale,
            xspec_settings=xspec_settings,
        )

    @classmethod
    def from_pyxspec_model(cls, model, spectrum_index=None):
        """
        Create a spectrum from a PyXspec model object.

        Parameters
        ----------
        model : :class:`~xspec.Model`
            The PyXspec model object.
        spectrum_index : integer, optional
            The index of the spectrum in the model object to
            use for getting the energies and fluxes. Default: Will
            use the first valid value of the spectral index, starting
            from zero.
        """
        if spectrum_index is None:
            spectrum_index = 0
        try:
            ebins = model.energies(spectrum_index)
        except Exception as e:
            if spectrum_index == 0:
                spectrum_index = 1
                ebins = model.energies(spectrum_index)
            else:
                raise e
        flux = model.values(spectrum_index)
        return cls._from_ext_model(ebins, flux)


class CountRateSpectrum(BaseSpectrum):
    _units = "photon/(s*keV)"

    @property
    def rate(self):
        return self._spec

    @property
    def luminosity(self):
        return self._energy_spec

    @property
    def total_rate(self):
        return self._total_spec

    @property
    def total_luminosity(self):
        return self._total_energy_spec

    @property
    def binned_rate(self):
        return self._binned_spec

    @property
    def binned_luminosity(self):
        return self._binned_energy_spec

    @property
    def rate_per_wavelength(self):
        return self.binned_rate / self.dwv

    @property
    def luminosity_per_wavelength(self):
        return self.binned_luminosity / self.dwv

    @property
    def rate_per_frequency(self):
        return self.binned_rate / self.df

    @property
    def luminosity_per_frequency(self):
        return self.binned_luminosity / self.df

    def generate_energies(self, t_exp, prng=None, quiet=False):
        """
        Generate photon energies from this count rate spectrum given an
        exposure time.

        Parameters
        ----------
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only
            be specified if you have a reason to generate the same
            set of random numbers, such as for a test. Default is None,
            which sets the seed based on the system time.
        quiet : boolean, optional
            If True, log messages will not be displayed when
            creating energies. Useful if you have to loop over
            a lot of spectra. Default: False
        """
        t_exp = parse_value(t_exp, "s")
        prng = parse_prng(prng)
        rate = self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng, self.binscale, quiet=quiet)
        energies = u.Quantity(energy, "keV")
        return energies

    @classmethod
    def from_spectrum(cls, spec, redshift=0.0, dist=None, cosmology=None):
        """
        Create a CountRateSpectrum from a Spectrum object, assuming that the latter
        is in the observer frame. To get the CountRateSpectrum in the source
        frame, you must either provide a redshift or a distance to the source to
        convert to photon count rate and luminosity. Note that the energy bins
        of the original Spectrum object will be cosmologically blueshifted if
        the redshift is nonzero.

        Parameters
        ----------
        spec : :class:`~soxs.spectra.base.Spectrum`
            The Spectrum object to create the CountRateSpectrum from.
        redshift : float
            The redshift to the source, assuming it is cosmological.
        dist : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The distance to the source, if it is not cosmological. If a unit
            is not specified, it is assumed to be in kpc.
        cosmology : :class:`~astropy.cosmology.Cosmology` object
            An AstroPy cosmology object used to determine the luminosity
            distance if needed. If not set, the default is the Planck 2018
            cosmology.
        """
        D_A = distance_handler(redshift, dist, cosmology)
        rate = spec.flux.value * (4.0 * np.pi * D_A**2 * (1.0 + redshift) ** 2)
        return cls(spec.ebins.value * (1.0 + redshift), rate, binscale=spec.binscale)

    @classmethod
    def from_powerlaw(cls, photon_index, norm, emin, emax, nbins, binscale="linear"):
        """
        Create a count rate spectrum from a power-law model, in the frame of the
        source. The form of the model is R(E) = norm*(e/(1 keV))**-photon_index.

        Parameters
        ----------
        photon_index : float
            The photon index of the source.
        norm : float
            The normalization of the source in units of
            the spectrum, photons/s/keV, at 1 keV
            in the source frame.
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy of the spectrum in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy of the spectrum in keV.
        nbins : integer
            The number of bins in the spectrum.
        binscale : string, optional
            The scale of the energy binning: "linear" or "log".
            Default: "linear"
        """
        return super().from_powerlaw(photon_index, 0.0, norm, emin, emax, nbins, binscale=binscale)

    @classmethod
    def from_xstar_model(cls, xstar_file, which_spectrum):
        """
        Load a count rate spectrum from an XSTAR output FITS file.
        See the XSTAR manual (link below) for details on how to generate the
        output file.

        https://heasarc.gsfc.nasa.gov/docs/software/xstar/

        Parameters
        ----------
        xstar_file : string
            The path to the XSTAR output FITS file, typically with
            "spect1" in the filename.
        which_spectrum : string
            Which spectrum to extract from the file. Options are:
            "emit_outward", "emit_inward", "transmitted", "incident"
        """
        from astropy.io import fits

        with fits.open(xstar_file) as f:
            if "XSTAR_SPECTRA" not in f:
                raise ValueError("No XSTAR_SPECTRA extension found in file!")
            hdu = f["XSTAR_SPECTRA"]
            # energies in file are bin left edges in eV
            elow = hdu.data["energy"].astype("float64") * 1.0e-3  # convert from eV to keV
            # bin sizes are logarithmic, use the size of the last bin to get the upper edge
            de = elow[-1] / elow[-2]
            emax = elow[-1] * de
            ebins = np.append(elow, emax)
            emid = 0.5 * (ebins[1:] + ebins[:-1])
            if which_spectrum not in hdu.data.dtype.names:
                raise ValueError(
                    f"Spectrum '{which_spectrum}' not found in XSTAR file! "
                    f"Options are: {hdu.data.dtype.names}"
                )
            rate = hdu.data[which_spectrum] * 1.0e38
            rate /= emid * erg_per_keV
            rate *= erg_per_keV
            return cls(ebins, rate, binscale="custom")

    def get_lum_in_band(self, emin, emax):
        """
        Determine the total photon count rate and luminosity
        within a band specified by an energy range.

        Parameters
        ----------
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy in the band, in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy in the band, in keV.

        Returns
        -------
        A tuple of values for the rate/luminosity in the
        band: the first value is in terms of the photon
        rate, the second value is in terms of the energy rate
        (or simply the luminosity).
        """
        return self._get_total_in_band(emin, emax)

    def get_flux_in_band(self, emin, emax):
        issue_deprecation_warning(
            "'get_flux_in_band' is deprecated for CountRateSpectrum, use 'get_lum_in_band' instead!",
        )
        return self.get_lum_in_band(emin, emax)


class ConvolvedSpectrum(CountRateSpectrum):
    def __init__(self, ebins, flux, arf, rmf=None, binscale="linear", noisy=False, exp_time=None):
        from numbers import Number

        from soxs.response import (
            AuxiliaryResponseFile,
            FlatResponse,
            RedistributionMatrixFile,
        )

        super().__init__(ebins, flux, binscale=binscale)
        if isinstance(arf, Number):
            arf = FlatResponse(ebins[0], ebins[-1], arf, ebins.size - 1)
        elif isinstance(arf, str):
            arf = AuxiliaryResponseFile(arf)
        self.arf = arf
        if rmf is not None and isinstance(rmf, str):
            rmf = RedistributionMatrixFile(rmf)
        self.rmf = rmf
        self.noisy = noisy
        if exp_time is not None:
            exp_time = parse_value(exp_time, "s") * u.s
        self.exp_time = exp_time
        if self.exp_time is not None and self.noisy:
            self._spec_err = (
                np.sqrt(self._spec * self.exp_time * self.de).value * u.photon / (self.exp_time * self.de)
            )
        else:
            self._spec_err = None

    _counts = None

    @property
    def counts(self):
        if self._counts is None and self.exp_time is not None:
            counts = (self.rate * self.exp_time).value
            if self.noisy:
                counts = np.rint(counts.value).astype("int")
            self._counts = counts * u.photon
        return self._counts

    _counts_err = None

    @property
    def counts_err(self):
        if self._spec_err is not None:
            self._counts_err = self._spec_err * self.exp_time * self.de
        elif self._counts_err is None and self.exp_time is not None:
            self._counts_err = np.sqrt(self.counts.value) * u.photon
        return self._counts_err

    def _check_binning_units(self, other):
        super()._check_binning_units(other)
        if self.noisy != other.noisy:
            raise RuntimeError("The noisy flags for these two spectra are not the same!")
        bad_exp = self.exp_time is None and other.exp_time is not None
        bad_exp |= self.exp_time is not None and other.exp_time is None
        if self.exp_time is not None and other.exp_time is not None:
            bad_exp |= not np.isclose(self.exp_time.value, other.exp_time.value)
        if bad_exp:
            raise RuntimeError("The exposure times for these two spectra are not the same!")
        bad_rmf = self.rmf is not None and other.rmf is None
        bad_rmf |= self.rmf is None and other.rmf is not None
        bad_rmf |= self.rmf.filename != other.rmf.filename
        if bad_rmf:
            raise RuntimeError("The RMF for these two spectra are not the same!")
        bad_arf = self.arf is not None and other.arf is None
        bad_arf |= self.arf is None and other.arf is not None
        bad_arf |= self.arf.filename != other.arf.filename
        if bad_arf:
            raise RuntimeError("The ARF for these two spectra are not the same!")

    def __copy__(self):
        # Create a new instance and copy attributes
        return type(self)(
            self.ebins,
            self._spec,
            self.arf,
            rmf=self.rmf,
            binscale=self.binscale,
            noisy=self.noisy,
            exp_time=self.exp_time,
        )

    def __deepcopy__(self, memo):
        # Create a new instance and deep copy attributes
        from soxs.response import (
            AuxiliaryResponseFile,
            FlatResponse,
            RedistributionMatrixFile,
        )

        if self.arf is None:
            arf = None
        elif isinstance(self.arf, FlatResponse):
            arf = FlatResponse(self.arf.ebins[0], self.arf.ebins[-1], self.arf.max_area, self.arf.nbins)
        else:
            arf = AuxiliaryResponseFile(self.arf.filename)
        rmf = None if self.rmf is None else RedistributionMatrixFile(self.rmf.filename)
        return type(self)(
            copy.deepcopy(self.ebins, memo),
            copy.deepcopy(self._spec, memo),
            arf,
            rmf=rmf,
            binscale=self.binscale,
            noisy=self.noisy,
            exp_time=self.exp_time,
        )

    def __add__(self, other):
        self._check_binning_units(other)
        return ConvolvedSpectrum(
            self.ebins,
            self._spec + other._spec,
            self.arf,
            binscale=self.binscale,
            rmf=self.rmf,
            noisy=self.noisy,
            exp_time=self.exp_time,
        )

    def __sub__(self, other):
        self._check_binning_units(other)
        return ConvolvedSpectrum(
            self.ebins,
            self._spec - other._spec,
            self.arf,
            binscale=self.binscale,
            rmf=self.rmf,
            noisy=self.noisy,
            exp_time=self.exp_time,
        )

    def __mul__(self, other):
        return ConvolvedSpectrum(
            self.ebins,
            other * self._spec,
            self.arf,
            binscale=self.binscale,
            rmf=self.rmf,
            noisy=self.noisy,
            exp_time=self.exp_time,
        )

    @classmethod
    def from_pha_file(cls, specfile):
        """
        Create a :class:`~soxs.spectra.base.ConvolvedSpectrum` object by reading it in
        from a PHA/PI file.

        Parameters
        ----------
        specfile : string
            The PHA/PI file to be opened.
        """
        from astropy.io import fits

        from soxs.response import AuxiliaryResponseFile, RedistributionMatrixFile

        with fits.open(specfile) as f:
            hdu = f["SPECTRUM"]
            if "COUNT_RATE" in hdu.data:
                rate = hdu.data["COUNT_RATE"].astype("float64")
            elif "RATE" in hdu.data:
                rate = hdu.data["RATE"].astype("float64")
            elif "COUNTS" in hdu.data and "EXPOSURE" in hdu.header:
                rate = hdu.data["COUNTS"].astype("float64") / hdu.header["EXPOSURE"]
            else:
                raise RuntimeError("Cannot determine count rate from this spectrum file!!")
            rate *= u.photon / u.s
            arf = AuxiliaryResponseFile(hdu.header["ANCRFILE"])
            rmf = RedistributionMatrixFile(hdu.header["RESPFILE"])
        ebins = np.append(rmf.ebounds_data["E_MIN"], rmf.ebounds_data["E_MAX"][-1]) * u.keV
        de = np.diff(ebins)
        rate /= de
        return cls(
            ebins,
            rate,
            arf,
            binscale="linear",
            rmf=rmf,
            noisy=bool(hdu.header["POISSERR"]),
            exp_time=hdu.header["EXPOSURE"],
        )

    def to_pha_file(self, specfile, overwrite=False):
        """
        Write a :class:`~soxs.spectra.base.ConvolvedSpectrum` object to a
        PHA/PI file.

        Parameters
        ----------
        specfile : string
            The PHA/PI file to be written.
        overwrite : boolean, optional
            Whether to overwrite an already existing file. Default: False
        """
        from soxs.events.spectra import _write_spectrum

        event_params = {
            "RESPFILE": os.path.split(self.rmf.filename)[-1],
            "ANCRFILE": os.path.split(self.arf.filename)[-1],
            "TELESCOP": self.rmf.header["TELESCOP"],
            "INSTRUME": self.rmf.header["INSTRUME"],
            "MISSION": self.rmf.header.get("MISSION", ""),
            "CHANTYPE": self.rmf.chan_type,
        }
        if self.exp_time is not None:
            event_params["EXPOSURE"] = self.exp_time.value
            spec = self.counts.value
        else:
            spec = self.rate.value

        bins = (np.arange(self.rmf.n_ch) + self.rmf.cmin).astype("int32")
        _write_spectrum(
            bins,
            spec,
            event_params,
            specfile,
            overwrite=overwrite,
            noisy=self.noisy,
        )

    @classmethod
    def convolve(cls, spectrum, arf, use_arf_energies=False, rmf=None):
        """
        Generate a model convolved spectrum by convolving a spectrum with an
        ARF, and optionally an RMF.

        Parameters
        ----------
        spectrum : :class:`~soxs.spectra.base.Spectrum` object
            The input spectrum to convolve with.
        arf : string or :class:`~soxs.response.AuxiliaryResponseFile`
            The ARF to use in the convolution.
        use_arf_energies : boolean, optional
            If True, use the energy binning of the ARF for
            the convolved spectrum. Default: False, which uses
            the original spectral binning.
        rmf : string or :class:`~soxs.response.RedistributionMatrixFile`, optional
            The RMF to use in the convolution if desired. Default is no RMF.
        """
        from soxs.response import AuxiliaryResponseFile, RedistributionMatrixFile

        if not isinstance(arf, AuxiliaryResponseFile):
            arf = AuxiliaryResponseFile(arf)
        if rmf is not None and not isinstance(rmf, RedistributionMatrixFile):
            rmf = RedistributionMatrixFile(rmf)
        if use_arf_energies or rmf is not None:
            rate = (
                regrid_spectrum(
                    arf.ebins,
                    spectrum.ebins.value,
                    spectrum._spec.value * spectrum.de.value,
                )
                * arf.eff_area
            )
            if rmf is not None:
                rate = rmf.convolve_spectrum((rate * arf.de).value, 1.0, noisy=False, rate=True)
            rate = u.Quantity(rate / arf.de, "keV-1 ph s-1")
            binscale = "linear"
            ebins = u.Quantity(arf.ebins, "keV")
        else:
            earea = arf.interpolate_area(spectrum.emid.value)
            rate = spectrum._spec * earea
            binscale = spectrum.binscale
            ebins = spectrum.ebins
        return cls(ebins, rate, arf, binscale=binscale, rmf=rmf)

    def new_spec_from_band(self, emin, emax):
        """
        Create a new :class:`~soxs.spectra.base.ConvolvedSpectrum` object
        from a subset of an existing one defined by a particular
        energy band.

        Parameters
        ----------
        emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The minimum energy of the band in keV.
        emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The maximum energy of the band in keV.
        """
        ebins, flux = self._new_spec_from_band(emin, emax)
        return ConvolvedSpectrum(
            ebins,
            flux,
            self.arf,
            binscale=self.binscale,
            noisy=self.noisy,
            rmf=self.rmf,
            exp_time=self.exp_time,
        )

    def deconvolve(self):
        """
        Return the deconvolved :class:`~soxs.spectra.base.Spectrum`
        object associated with this convolved spectrum.
        """
        if self.rmf is not None:
            raise NotImplementedError("deconvolve is not implemented for ConvolvedSpectra with RMFs!")
        earea = self.arf.interpolate_area(self.emid)
        flux = self._spec / earea
        flux = np.nan_to_num(flux.value)
        return Spectrum(self.ebins.value, flux, binscale=self.binscale)

    def generate_energies(self, t_exp, prng=None, quiet=False):
        """
        Generate photon energies from this convolved spectrum given an
        exposure time.

        Parameters
        ----------
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only
            be specified if you have a reason to generate the same
            set of random numbers, such as for a test. Default is None,
            which sets the seed based on the system time.
        quiet : boolean, optional
            If True, log messages will not be displayed when
            creating energies. Useful if you have to loop over
            a lot of spectra. Default: False
        """
        if self.rmf is not None:
            raise NotImplementedError("generate_energies is not implemented for ConvolvedSpectra with RMFs!")
        t_exp = parse_value(t_exp, "s")
        prng = parse_prng(prng)
        rate = self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng, self.binscale, quiet=quiet)
        earea = self.arf.interpolate_area(energy).value
        flux = np.sum(energy) * erg_per_keV / t_exp / earea.sum()
        energies = Energies(energy, flux)
        return energies

    @classmethod
    def from_constant(cls, const_value, emin, emax, nbins, binscale="linear"):
        raise NotImplementedError

    @classmethod
    def from_powerlaw(cls, photon_index, norm, emin, emax, nbins, binscale="linear"):
        raise NotImplementedError

    @classmethod
    def from_xspec_model(
        cls,
        model_string,
        params,
        emin,
        emax,
        nbins,
        binscale="linear",
        xspec_settings=None,
    ):
        raise NotImplementedError

    @classmethod
    def from_xspec_script(cls, infile, emin, emax, nbins, binscale="linear", xspec_settings=None):
        raise NotImplementedError

    def apply_foreground_absorption(self, nH, model="wabs", redshift=0.0, abund_table="angr"):
        raise NotImplementedError

    def add_absorption_line(self, line_center, line_width, equiv_width, line_type="gaussian"):
        raise NotImplementedError

    def add_absorption_lines(self, line_center, line_width, equiv_lines):
        raise NotImplementedError
