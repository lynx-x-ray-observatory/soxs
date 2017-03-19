import numpy as np
from soxs.spectra import Spectrum, ConvolvedSpectrum, \
    _generate_energies, Energies
from soxs.constants import erg_per_keV
from soxs.utils import parse_prng

class BackgroundSpectrum(Spectrum):
    _units = "photon/(cm**2*s*keV*arcmin**2)"
    def __init__(self, filename):
        emid, flux = np.loadtxt(filename, unpack=True)
        de = np.diff(emid)[0]
        ebins = np.append(emid-0.5*de, emid[-1]+0.5*de)
        super(BackgroundSpectrum, self).__init__(ebins, flux)

    def generate_energies(self, t_exp, area, fov, prng=None):
        """
        Generate photon energies from this background 
        spectrum given an exposure time, effective area, 
        and field of view.

        Parameters
        ----------
        t_exp : float
            The exposure time in seconds.
        area : float
            The effective area in cm**2. If one is creating 
            events for a SIMPUT file, a constant should be 
            used and it must be large enough so that a 
            sufficiently large sample is drawn for the ARF.
        fov : float
            The width of the field of view on a side in 
            arcminutes.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only 
            be specified if you have a reason to generate the same 
            set of random numbers, such as for a test. Default is None, 
            which sets the seed based on the system time. 
        """
        prng = parse_prng(prng)
        rate = area*fov*fov*self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng)
        flux = np.sum(energy)*erg_per_keV/t_exp/area
        energies = Energies(energy, flux)
        return energies

class ConvolvedBackgroundSpectrum(ConvolvedSpectrum):
    _units = "photon/(s*keV*arcmin**2)"

    def generate_energies(self, t_exp, fov, prng=None):
        """
        Generate photon energies from this convolved 
        background spectrum given an exposure time and 
        field of view.

        Parameters
        ----------
        t_exp : float
            The exposure time in seconds.
        fov : float
            The width of the field of view on a side 
            in arcminutes.
        prng : :class:`~numpy.random.RandomState` object, integer, or None
            A pseudo-random number generator. Typically will only 
            be specified if you have a reason to generate the same 
            set of random numbers, such as for a test. Default is None, 
            which sets the seed based on the system time. 
        """
        prng = parse_prng(prng)
        rate = fov*fov*self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng)
        earea = self.arf.interpolate_area(energy).value
        flux = np.sum(energy)*erg_per_keV/t_exp/earea.sum()
        energies = Energies(energy, flux)
        return energies
