from __future__ import division

import numpy as np
import subprocess
import tempfile
import shutil
import os
from soxs.utils import soxs_files_path, mylog
from soxs.cutils import broaden_lines
from soxs.constants import erg_per_keV, hc, \
    cosmic_elem, metal_elem, atomic_weights, clight, \
    m_u
import astropy.io.fits as pyfits
import astropy.units as u

class Energies(u.Quantity):
    def __new__(cls, energy, flux, flux_units):
        ret = u.Quantity.__new__(cls, energy, unit="keV")
        ret.flux = u.Quantity(flux, flux_units)
        return ret

def _generate_energies(spec, t_exp, rate, prng=None):
    cumspec = spec.cumspec
    n_ph = np.modf(t_exp*rate)
    n_ph = np.int64(n_ph[1]) + np.int64(n_ph[0] >= prng.uniform())
    mylog.info("Creating %d events from this spectrum." % n_ph)
    randvec = prng.uniform(size=n_ph)
    randvec.sort()
    return np.interp(randvec, cumspec, spec.ebins.value)

class Spectrum(object):
    _units = "photon/(cm**2*s*keV)"
    def __init__(self, ebins, flux):
        self.ebins = u.Quantity(ebins, "keV")
        self.emid = 0.5*(self.ebins[1:]+self.ebins[:-1])
        self.flux = u.Quantity(flux, self._units)
        self.nbins = len(self.emid)
        self.de = self.ebins[1]-self.ebins[0]
        self._compute_total_flux()

    def _compute_total_flux(self):
        self.total_flux = self.flux.sum()*self.de
        self.total_energy_flux = (self.flux*self.emid.to("erg")).sum()*self.de/(1.0*u.photon)
        cumspec = np.cumsum(self.flux.value*self.de.value)
        cumspec = np.insert(cumspec, 0, 0.0)
        cumspec /= cumspec[-1]
        self.cumspec = cumspec

    def __add__(self, other):
        if self.nbins != other.nbins or \
            not np.isclose(self.ebins.value, other.ebins.value).all():
            raise RuntimeError("Energy binning for these two "
                               "spectra is not the same!!")
        return Spectrum(self.ebins, self.flux+other.flux)

    def __mul__(self, other):
        return Spectrum(self.ebins, other*self.flux)

    __rmul__ = __mul__

    def __div__(self, other):
        return Spectrum(self.ebins, self.flux/other)

    def __truediv__(self, other):
        return Spectrum(self.ebins, self.flux/other)

    def __repr__(self):
        s = "Spectrum (%s - %s)\n" % (self.ebins[0], self.ebins[-1])
        s += "    Total Flux:\n    %s\n    %s\n" % (self.total_flux, self.total_energy_flux)
        return s

    def get_flux_in_band(self, emin, emax):
        """
        Determine the total flux within a band specified by an energy range. 

        Parameters
        ----------
        emin : float
            The minimum energy in the band, in keV.
        emax : float
            The maximum energy in the band, in keV.

        Returns
        -------
        A tuple of values for the flux/intensity in the band: the first value is in
        terms of the photon rate, the second value is in terms of the energy rate. 
        """
        range = np.logical_and(self.emid.value >= emin, self.emid.value <= emax)
        pflux = self.flux[range].sum()*self.de
        eflux = (self.flux*self.emid.to("erg"))[range].sum()*self.de/(1.0*u.photon)
        return pflux, eflux

    @classmethod
    def from_xspec_script(cls, infile, emin=0.01, emax=50.0, nbins=10000):
        """
        Create a model spectrum using a script file as input to XSPEC.

        Parameters
        ----------
        infile : string
            Path to the script file to use. 
        emin : float, optional
            The minimum energy of the spectrum in keV. Default: 0.01
        emax : float, optional
            The maximum energy of the spectrum in keV. Default: 50.0
        nbins : integer, optional
            The number of bins in the spectrum. Default: 10000
        """
        f = open(infile, "r")
        xspec_in = f.readlines()
        f.close()
        return cls._from_xspec(xspec_in, emin, emax, nbins)

    @classmethod
    def from_xspec_model(cls, model_string, params, emin=0.01, emax=50.0,
                         nbins=10000):
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
        emin : float, optional
            The minimum energy of the spectrum in keV. Default: 0.01
        emax : float, optional
            The maximum energy of the spectrum in keV. Default: 50.0
        nbins : integer, optional
            The number of bins in the spectrum. Default: 10000
        """
        xspec_in = []
        model_str = "%s &" % model_string
        for param in params:
            model_str += " %g &" % param
        model_str += " /*"
        xspec_in.append("model %s\n" % model_str)
        return cls._from_xspec(xspec_in, emin, emax, nbins)

    @classmethod
    def from_xspec(cls, model_string, params, emin=0.01, emax=50.0,
                   nbins=10000):
        mylog.warning("The 'from_xspec' method has been deprecated: "
                      "use 'from_xspec_model' instead.")
        cls.from_xspec_model(model_string, params, emin=emin, emax=emax, nbins=nbins)

    @classmethod
    def _from_xspec(cls, xspec_in, emin, emax, nbins):
        tmpdir = tempfile.mkdtemp()
        curdir = os.getcwd()
        os.chdir(tmpdir)
        xspec_in.append("dummyrsp %g %g %d lin\n" % (emin, emax, nbins))
        xspec_in += ["set fp [open spec_therm.xspec w+]\n",
                     "tclout energies\n", "puts $fp $xspec_tclout\n",
                     "tclout modval\n", "puts $fp $xspec_tclout\n",
                     "close $fp\n", "quit\n"]
        f_xin = open("xspec.in", "w")
        f_xin.writelines(xspec_in)
        f_xin.close()
        logfile = os.path.join(curdir, "xspec.log")
        with open(logfile, "ab") as xsout:
            subprocess.call(["xspec", "-", "xspec.in"],
                            stdout=xsout, stderr=xsout)
        f_s = open("spec_therm.xspec", "r")
        lines = f_s.readlines()
        f_s.close()
        ebins = np.array(lines[0].split()).astype("float64")
        de = np.diff(ebins)[0]
        flux = np.array(lines[1].split()).astype("float64")/de
        os.chdir(curdir)
        shutil.rmtree(tmpdir)
        return cls(ebins, flux)

    @classmethod
    def from_powerlaw(cls, photon_index, redshift, norm,
                      emin=0.01, emax=50.0, nbins=10000):
        """
        Create a spectrum from a power-law model.

        Parameters
        ----------
        photon_index : float
            The photon index of the source.
        redshift : float
            The redshift of the source.
        norm : float
            The normalization of the source in units of
            photons/s/cm**2/keV at 1 keV in the source frame.
        emin : float, optional
            The minimum energy of the spectrum in keV. Default: 0.01
        emax : float, optional
            The maximum energy of the spectrum in keV. Default: 50.0
        nbins : integer, optional
            The number of bins in the spectrum. Default: 10000
        """
        ebins = np.linspace(emin, emax, nbins+1)
        emid = 0.5*(ebins[1:]+ebins[:-1])
        flux = norm*(emid*(1.0+redshift))**(-photon_index)
        return cls(ebins, flux)

    @classmethod
    def from_file(cls, filename):
        """
        Read a spectrum from an ASCII text file. Accepts a file
        with two columns, the first being the center energy of the
        bin in keV and the second being the flux in photons/s/cm**2/keV, 
        assuming a linear binning with constant bin widths:

        Parameters
        ----------
        filename : string
            The path to the file containing the spectrum.
        """
        emid, flux = np.loadtxt(filename, unpack=True)
        de = np.diff(emid)[0]
        ebins = np.append(emid-0.5*de, emid[-1]+0.5*de)
        return cls(ebins, flux)

    @classmethod
    def from_constant(cls, const_flux, emin=0.01, emax=50.0, nbins=10000):
        """
        Create a spectrum from a constant model using XSPEC.

        Parameters
        ----------
        const_flux : float
            The value of the constant flux in the units of the spectrum. 
        emin : float, optional
            The minimum energy of the spectrum in keV. Default: 0.01
        emax : float, optional
            The maximum energy of the spectrum in keV. Default: 50.0
        nbins : integer, optional
            The number of bins in the spectrum. Default: 10000
        """
        ebins = np.linspace(emin, emax, nbins+1)
        flux = const_flux*np.ones(nbins)
        return cls(ebins, flux)

    def rescale_flux(self, new_flux, emin=None, emax=None, flux_type="photons"):
        """
        Rescale the flux of the spectrum, optionally using a
        specific energy band.

        Parameters
        ----------
        new_flux : float
            The new flux in units of photons/s/cm**2.
        emin : float, optional
            The minimum energy of the band to consider, in keV.
            Default: Use the minimum energy of the entire spectrum.
        emax : float, optional
            The maximum energy of the band to consider, in keV.
            Default: Use the maximum energy of the entire spectrum.
        flux_type : string, optional
            The units of the flux to use in the rescaling:
                "photons": photons/s/cm**2
                "energy": erg/s/cm**2
        """
        if emin is None:
            emin = self.ebins[0].value
        if emax is None:
            emax = self.ebins[-1].value
        idxs = np.logical_and(self.emid.value >= emin, self.emid.value <= emax)
        if flux_type == "photons":
            f = self.flux[idxs].sum()*self.de
        elif flux_type == "energy":
            f = (self.flux*self.emid.to("erg"))[idxs].sum()*self.de
        self.flux *= new_flux/f.value
        self._compute_total_flux()

    def write_file(self, specfile, clobber=False):
        """
        Write the spectrum to a file.

        Parameters
        ----------
        specfile : string
            The filename to write the file to.
        clobber : boolean, optional
            Whether or not to overwrite an existing file. Default: False
        """
        if os.path.exists(specfile) and not clobber:
            raise IOError("File %s exists and clobber=False!" % specfile)
        header = "Energy\tFlux\nkeV\t%s" % self._units
        np.savetxt(specfile, np.transpose([self.emid, self.flux]), 
                   delimiter="\t", header=header)

    def apply_foreground_absorption(self, nH):
        """
        Given a hydrogen column density, apply
        galactic foreground absorption to the spectrum. 

        Parameters
        ----------
        nH : float
            The hydrogen column in units of 10**22 atoms/cm**2
        """
        sigma = wabs_cross_section(self.emid.value)
        self.flux *= np.exp(-nH*1.0e22*sigma)
        self._compute_total_flux()

    def generate_energies(self, t_exp, area, prng=None):
        """
        Generate photon energies from this spectrum given an
        exposure time and effective area.

        Parameters
        ----------
        t_exp : float
            The exposure time in seconds.
        area : float
            The effective area in cm**2, If one is creating events for a SIMPUT file, 
            a constant should be used and it must be large enough  so that a sufficiently 
            large sample is drawn for the ARF.
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        """
        if prng is None:
            prng = np.random
        rate = area*self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng=prng)
        flux = np.sum(energy)*erg_per_keV/t_exp/area
        energies = Energies(energy, flux, "erg/(cm**2*s)")
        return energies

class ApecGenerator(object):
    r"""
    Initialize a thermal gas emission model from the AtomDB APEC tables
    available at http://www.atomdb.org. This code borrows heavily from Python
    routines used to read the APEC tables developed by Adam Foster at the
    CfA (afoster@cfa.harvard.edu).

    Parameters
    ----------
    emin : float
        The minimum energy for the spectral model.
    emax : float
        The maximum energy for the spectral model.
    nbins : integer
        The number of bins in the spectral model.
    apec_root : string, optional
        The directory root where the APEC model files are stored. If 
        not provided, the default is to grab them from the tables stored
        with SOX.
    apec_vers : string, optional
        The version identifier string for the APEC files, e.g.
        "2.0.2"
    broadening : boolean, optional
        Whether or not the spectral lines should be thermally
        and velocity broadened.

    Examples
    --------
    >>> apec_model = ApecGenerator(0.05, 50.0, 1000, apec_vers="3.0.3",
    ...                            broadening=True)
    """
    def __init__(self, emin, emax, nbins, apec_root=None,
                 apec_vers="2.0.2", broadening=True):
        self.emin = emin
        self.emax = emax
        self.nbins = nbins
        self.ebins = np.linspace(self.emin, self.emax, nbins+1)
        self.de = np.diff(self.ebins)
        self.emid = 0.5*(self.ebins[1:]+self.ebins[:-1])
        if apec_root is None:
            apec_root = soxs_files_path
        self.cocofile = os.path.join(apec_root, "apec_v%s_coco.fits" % apec_vers)
        self.linefile = os.path.join(apec_root, "apec_v%s_line.fits" % apec_vers)
        if not os.path.exists(self.cocofile) or not os.path.exists(self.linefile):
            raise IOError("Cannot find the APEC files!\n %s\n, %s" % (self.cocofile,
                                                                      self.linefile))
        self.wvbins = hc/self.ebins[::-1]
        self.broadening = broadening
        try:
            self.line_handle = pyfits.open(self.linefile)
        except IOError:
            raise IOError("LINE file %s does not exist" % self.linefile)
        try:
            self.coco_handle = pyfits.open(self.cocofile)
        except IOError:
            raise IOError("COCO file %s does not exist" % self.cocofile)

        self.Tvals = self.line_handle[1].data.field("kT")
        self.nT = len(self.Tvals)
        self.dTvals = np.diff(self.Tvals)
        self.minlam = self.wvbins.min()
        self.maxlam = self.wvbins.max()

    def _make_spectrum(self, kT, element, velocity, line_fields, 
                       coco_fields, scale_factor):

        tmpspec = np.zeros(self.nbins)

        i = np.where((line_fields['element'] == element) &
                     (line_fields['lambda'] > self.minlam) &
                     (line_fields['lambda'] < self.maxlam))[0]

        E0 = hc/line_fields['lambda'][i].astype("float64")*scale_factor
        amp = line_fields['epsilon'][i].astype("float64")
        if self.broadening:
            sigma = 2.*kT*erg_per_keV/(atomic_weights[element]*m_u)
            sigma += 2.0*velocity*velocity
            sigma = E0*np.sqrt(sigma)/clight
            vec = broaden_lines(E0, sigma, amp, self.ebins)
        else:
            vec = np.histogram(E0, self.ebins, weights=amp)[0]
        tmpspec += vec

        ind = np.where((coco_fields['Z'] == element) &
                       (coco_fields['rmJ'] == 0))[0]
        if len(ind) == 0:
            return tmpspec
        else:
            ind = ind[0]

        n_cont = coco_fields['N_Cont'][ind]
        e_cont = coco_fields['E_Cont'][ind][:n_cont]*scale_factor
        continuum = coco_fields['Continuum'][ind][:n_cont]

        tmpspec += np.interp(self.emid, e_cont, continuum)*self.de/scale_factor

        n_pseudo = coco_fields['N_Pseudo'][ind]
        e_pseudo = coco_fields['E_Pseudo'][ind][:n_pseudo]*scale_factor
        pseudo = coco_fields['Pseudo'][ind][:n_pseudo]

        tmpspec += np.interp(self.emid, e_pseudo, pseudo)*self.de/scale_factor

        return tmpspec*scale_factor

    def _preload_data(self, index):
        line_data = self.line_handle[index+2].data
        coco_data = self.coco_handle[index+2].data
        line_fields = ('element', 'lambda', 'epsilon')
        coco_fields = ('Z', 'rmJ', 'N_Cont', 'E_Cont', 'Continuum', 'N_Pseudo',
                       'E_Pseudo', 'Pseudo')
        line_fields = {el: line_data.field(el) for el in line_fields}
        coco_fields = {el: coco_data.field(el) for el in coco_fields}
        return line_fields, coco_fields

    def get_spectrum(self, kT, abund, redshift, norm, velocity=0.0):
        """
        Get a thermal emission spectrum.

        Parameters
        ----------
        kT : float
            The temperature in keV.
        abund : float
            The metal abundance in solar units. 
        redshift : float
            The redshift.
        norm : float
            The normalization of the model, in the standard Xspec units of
            1.0e-14*EM/(4*pi*(1+z)**2*D_A**2).
        velocity : float, optional
            The velocity broadening parameter, in units of km/s. Default: 0.0
        """
        v = velocity*1.0e5
        tindex = np.searchsorted(self.Tvals, kT)-1
        if tindex >= self.Tvals.shape[0]-1 or tindex < 0:
            return np.zeros(self.nbins)
        cspec = np.zeros((2, self.nbins))
        mspec = np.zeros((2, self.nbins))
        scale_factor = 1./(1.+redshift)
        dT = (kT-self.Tvals[tindex])/self.dTvals[tindex]
        for i, ikT in enumerate([tindex, tindex+1]):
            line_fields, coco_fields = self._preload_data(ikT)
            # First do H,He, and trace elements
            for elem in cosmic_elem:
                cspec[i,:] += self._make_spectrum(self.Tvals[ikT], elem, v, line_fields, 
                                                  coco_fields, scale_factor)
            # Next do the metals
            for elem in metal_elem:
                mspec[i,:] += self._make_spectrum(self.Tvals[ikT], elem, v, line_fields, 
                                                  coco_fields, scale_factor)
        cosmic_spec = cspec[0,:]*(1.-dT)+cspec[1,:]*dT
        metal_spec = mspec[0,:]*(1.-dT)+mspec[1,:]*dT
        spec = 1.0e14*norm*(cosmic_spec + abund*metal_spec)/self.de
        return Spectrum(self.ebins, spec)

def wabs_cross_section(E):
    emax = np.array([0.0, 0.1, 0.284, 0.4, 0.532, 0.707, 0.867, 1.303, 1.840, 
                     2.471, 3.210, 4.038, 7.111, 8.331, 10.0])
    c0 = np.array([17.3, 34.6, 78.1, 71.4, 95.5, 308.9, 120.6, 141.3,
                   202.7,342.7,352.2,433.9,629.0,701.2])
    c1 = np.array([608.1, 267.9, 18.8, 66.8, 145.8, -380.6, 169.3,
                   146.8, 104.7, 18.7, 18.7, -2.4, 30.9, 25.2]) 
    c2 = np.array([-2150., -476.1 ,4.3, -51.4, -61.1, 294.0, -47.7,
                   -31.5, -17.0, 0.0, 0.0, 0.75, 0.0, 0.0])
    idxs = np.minimum(np.searchsorted(emax, E)-1, 13)
    sigma = (c0[idxs]+c1[idxs]*E+c2[idxs]*E*E)*1.0e-24/E**3
    return sigma

def get_wabs_absorb(e, nH):
    sigma = wabs_cross_section(e)
    return np.exp(-nH*1.0e22*sigma)

class ConvolvedSpectrum(Spectrum):
    _units = "photon/(s*keV)"
    def __init__(self, spectrum, arf):
        """
        Generate a convolved spectrum by convolving a spectrum with an
        ARF.

        Parameters
        ----------
        spectrum : :class:`~soxs.spectra.Spectrum` object
            The input spectrum to convolve with.
        arf : string or :class:`~soxs.instrument.AuxiliaryResponseFile`
            The ARF to use in the convolution.
        """
        from soxs.instrument import AuxiliaryResponseFile
        if not isinstance(arf, AuxiliaryResponseFile):
            arf = AuxiliaryResponseFile(arf)
        self.arf = arf
        earea = arf.interpolate_area(spectrum.emid)
        rate = spectrum.flux * earea
        super(ConvolvedSpectrum, self).__init__(spectrum.ebins, rate)

    def generate_energies(self, t_exp, prng=None):
        """
        Generate photon energies from this convolved spectrum given an
        exposure time.

        Parameters
        ----------
        t_exp : float
            The exposure time in seconds.
        prng : :class:`~numpy.random.RandomState` object or :mod:`~numpy.random`, optional
            A pseudo-random number generator. Typically will only be specified
            if you have a reason to generate the same set of random numbers, such as for a
            test. Default is the :mod:`numpy.random` module.
        """
        if prng is None:
            prng = np.random
        rate = self.total_flux.value
        energy = _generate_energies(self, t_exp, rate, prng=prng)
        flux = np.sum(energy)*erg_per_keV/t_exp
        energies = Energies(energy, flux, "erg/s")
        return energies
