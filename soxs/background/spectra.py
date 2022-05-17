import numpy as np
from soxs.spectra import Spectrum, ConvolvedSpectrum, \
    _generate_energies, Energies
from soxs.constants import erg_per_keV
from soxs.utils import parse_prng, parse_value
from soxs.response import AuxiliaryResponseFile


class BackgroundSpectrum(Spectrum):
    _units = "photon/(cm**2*s*keV*arcmin**2)"

    def __init__(self, ebins, flux, binscale="linear"):
        super().__init__(ebins, flux, binscale=binscale)

    @classmethod
    def from_spectrum(cls, spec, fov):
        """
        Create a background spectrum from a regular
        :class:`~soxs.spectra.Spectrum` object and the width
        of a field of view on a side.

        Parameters
        ----------
        spec : :class:`~soxs.spectra.Spectrum`
            The spectrum to be used.
        fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The width of the field of view on a side in 
            arcminutes.
        """
        fov = parse_value(fov, "arcmin")
        flux = spec.flux.value/fov/fov
        return cls(spec.ebins.value, flux, binscale=spec.binscale)

    def generate_energies(self, t_exp, area, fov, prng=None, 
                          quiet=False):
        """
        Generate photon energies from this background 
        spectrum given an exposure time, effective area, 
        and field of view.

        Parameters
        ----------
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        area : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The effective area in cm**2. If one is creating 
            events for a SIMPUT file, a constant should be 
            used and it must be large enough so that a 
            sufficiently large sample is drawn for the ARF.
        fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The width of the field of view on a side in 
            arcminutes.
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
        fov = parse_value(fov, "arcmin")
        area = parse_value(area, "cm**2")
        prng = parse_prng(prng)
        rate = area*fov*fov*self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng, self.binscale, quiet=quiet)
        flux = np.sum(energy)*erg_per_keV/t_exp/area
        energies = Energies(energy, flux)
        return energies

    def to_spectrum(self, fov):
        fov = parse_value(fov, "arcmin")
        flux = self.flux.value*fov*fov
        return Spectrum(self.ebins.value, flux, binscale=self.binscale)

    def __mul__(self, other):
        if isinstance(other, AuxiliaryResponseFile):
            return ConvolvedBackgroundSpectrum.convolve(self, other)
        else:
            return BackgroundSpectrum(self.ebins, other*self.flux, 
                                      binscale=self.binscale)

    __rmul__ = __mul__


class ConvolvedBackgroundSpectrum(ConvolvedSpectrum):
    _units = "photon/(s*keV*arcmin**2)"

    @classmethod
    def from_spectrum(cls, spec, fov):
        """
        Create a convolved background spectrum from a regular
        :class:`~soxs.spectra.ConvolvedSpectrum` object and the width
        of a field of view on a side.

        Parameters
        ----------
        spec : :class:`~soxs.spectra.ConvolvedSpectrum`
            The spectrum to be used.
        fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The width of the field of view on a side in 
            arcminutes.
        """
        fov = parse_value(fov, "arcmin")
        flux = spec.flux.value/fov/fov
        return cls(spec.ebins.value, flux, spec.binscale)

    def generate_energies(self, t_exp, fov, prng=None, 
                          quiet=False):
        """
        Generate photon energies from this convolved 
        background spectrum given an exposure time and 
        field of view.

        Parameters
        ----------
        t_exp : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The exposure time in seconds.
        fov : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The width of the field of view on a side 
            in arcminutes.
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
        fov = parse_value(fov, "arcmin")
        prng = parse_prng(prng)
        rate = fov*fov*self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng, self.binscale, quiet=quiet)
        earea = self.arf.interpolate_area(energy).value
        flux = np.sum(energy)*erg_per_keV/t_exp/earea.sum()
        energies = Energies(energy, flux)
        return energies
