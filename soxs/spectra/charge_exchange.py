import numpy as np
from tqdm.auto import tqdm

from soxs.constants import abund_tables, elem_names
from soxs.spectra import Spectrum
from soxs.utils import mylog, parse_value, soxs_cfg

elem_full = [
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    28,
]


class ACX2Generator:
    _one_ion = False
    """
    Generate charge exchange spectra using the ACX2 model. This
    class assumes that the balance of ions in the recombining
    plasma can be determined by an input temperature. To use this
    model, you must have the acx2 (https://acx2.readthedocs.io/)
    and pyatomdb (https://atomdb.readthedocs.io/) packages installed.

    Parameters
    ----------
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The minimum energy for the spectral model.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The maximum energy for the spectral model.
    nbins : integer
        The number of bins in the spectral model.
    collntype : integer, optional
        The type of collision parameter to use:
        1 - center of mass energy (kev/u)
        2 - center of mass velocity (km/s)
        3 - donor ion velocity (km/s)
        4 - recombining ion velocity (km/s)
        Default: 1
    acx_model : integer, optional
        ACX model to fall back on, from 1 to 8. Default: 8
    recomb_type : integer, optional
        The type of recombination to use:
        1 - single recombination
        2 - full recombination
        3 - full renormalized recombination
        Default: 1
    binscale : string, optional
        The scale of the energy binning: "linear" or "log".
        Default: "linear"
    var_elem : list of strings, optional
        The names of elements to allow to vary freely from the single
        abundance parameter. These must be strings like ["O", "N", "He"].
        Default: None
    abund_table : string or array_like, optional
        The abundance table to be used for solar abundances.
        Either a string corresponding to a built-in table or an array
        of 30 floats corresponding to the abundances of each element
        relative to the abundance of H. Default is set in the SOXS
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
        "cl17.03" : the abundance table used in Cloudy v17.03.
    """

    def __init__(
        self,
        emin,
        emax,
        nbins,
        collntype=1,
        acx_model=8,
        recomb_type=1,
        binscale="linear",
        var_elem=None,
        abund_table=None,
    ):
        try:
            from acx2 import ACXModel, __version__ as model_vers
        except ImportError:
            raise ImportError(
                "You must have the acx2 and pyatomdb packages "
                "installed to use the ACX2Generator class!"
            )
        self.model_vers = model_vers
        self.model = ACXModel()
        self.binscale = binscale
        emin = parse_value(emin, "keV")
        emax = parse_value(emax, "keV")
        self.emin = emin
        self.emax = emax
        self.nbins = nbins
        if binscale == "linear":
            self.ebins = np.linspace(self.emin, self.emax, nbins + 1)
        elif binscale == "log":
            self.ebins = np.logspace(
                np.log10(self.emin), np.log10(self.emax), nbins + 1
            )
        self.de = np.diff(self.ebins)
        self.emid = 0.5 * (self.ebins[1:] + self.ebins[:-1])
        self.elements = elem_full
        if var_elem is None:
            self.var_elem = np.empty(0, dtype="int")
        else:
            if len(var_elem) != len(set(var_elem)):
                raise RuntimeError(
                    'Duplicates were found in the "var_elem" list! %s' % var_elem
                )
            self.var_elem = [elem_names.index(e) for e in var_elem]
            self.var_elem.sort()
            self.var_elem = np.array(self.var_elem, dtype="int")
            self.var_elem_names = [elem_names[e] for e in self.var_elem]
            self.var_elem_idxs = np.array(
                [self.elements.index(e) for e in self.var_elem], dtype="int"
            )
        self.num_var_elem = len(self.var_elem)
        self.num_elements = len(self.elements)
        if abund_table is None:
            abund_table = soxs_cfg.get("soxs", "abund_table")
        if not isinstance(abund_table, str):
            if len(abund_table) != 30:
                raise RuntimeError(
                    "User-supplied abundance tables must be 30 elements long!"
                )
            self.atable = np.concatenate([[0.0], np.array(abund_table)])
        else:
            self.atable = abund_tables[abund_table].copy()
        self._atable = self.atable.copy()
        self._atable[1:] /= abund_tables["angr"][1:]
        self._atable = self._atable[self.elements]

        if len(self.model.DonorList) == 0:
            self.model.add_donor("H", elements=self.elements)
            self.model.add_donor("He", elements=self.elements)

        # acx fallback model
        self.acx_model = acx_model
        self.model.set_acxmodel(acx_model)

        # set recombination type (1 = single ['solar wind'], 2 = full ['comet'])
        self.recomb_type = recomb_type
        self.model.set_recombtype(recomb_type)

        # set collision type and units
        if collntype == 1:
            self.coll_units = "keV/u"
        else:
            self.coll_units = "km/s"

        self.collntype = collntype
        self.model.set_collisiontype(collntype, self.coll_units)

    def _get_spectrum(self, redshift, He_frac, collnpar, tbroad, velocity):

        # Shift the energy bins to the source frame
        self.model.set_ebins(self.ebins * (1.0 + redshift))

        # set the H & He fraction from HeFrac for the donors
        self.model.set_donorabund(["H", "He"], [1 - He_frac, He_frac])

        # collision parameter parsing
        collnpar = parse_value(collnpar, self.coll_units)

        # calculate the spectrum
        spec = self.model.calc_spectrum(collnpar, Tbroaden=tbroad, vbroaden=velocity)

        return spec

    def get_spectrum(
        self,
        kT,
        collnpar,
        abund,
        He_frac,
        redshift,
        norm,
        elem_abund=None,
        velocity=0.0,
        tbroad=0.0,
    ):
        """
        Get a charge exchange spectrum, given a temperature that sets the
        ionization balance of the recombining plasma.

        Parameters
        ----------
        kT : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The temperature in keV.
        collnpar : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The collision parameter. Units are determined by the value of
            the collntype parameter set in the class constructor.
        abund : float
            The metal abundance in solar units.
        He_frac : float
            Number fraction of donor which is He (remainder is H).
        redshift : float
            The redshift.
        norm : float
            The normalization of the model, in units of
            EM/(4*pi*(1+z)**2*D_A**2), where EM = int N_H^r N_(H+He)^d dV.
        elem_abund : dict of element name, float pairs, optional
            A dictionary of elemental abundances in solar
            units to vary freely of the abund parameter, e.g.
            {"O": 0.4, "N": 0.3, "He": 0.9}. Default: None
        velocity : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The velocity broadening parameter, in units of
            km/s. Default: 0.0
        tbroad : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The temperature broadening parameter, in units of keV. Default: 0.0
        """
        # check velocity and tbroad
        velocity = parse_value(velocity, "km/s")
        tbroad = parse_value(tbroad, "keV")

        # set the abundances
        if elem_abund is None:
            elem_abund = {}
        abunds = abund * np.ones(self.num_elements)
        if len(elem_abund) > 0:
            if set(elem_abund.keys()) != set(self.var_elem_names):
                msg = (
                    "The supplied set of abundances is not the same as "
                    "that which was originally set!\nFree elements: "
                    f"{set(elem_abund.keys())}\n"
                    f"Abundances: {set(self.var_elem_names)}"
                )
                raise RuntimeError(msg)
            for k, v in elem_abund.items():
                idx = self.var_elem_names.index(k)
                abunds[self.var_elem_idxs[idx]] = v
        self.model.set_abund(abunds * self._atable, elements=self.elements)

        # set the temperature
        self.model.set_temperature(kT)

        spec = self._get_spectrum(redshift, He_frac, collnpar, tbroad, velocity)

        return Spectrum(self.ebins, norm * spec / self.de, binscale=self.binscale)


class OneACX2Generator(ACX2Generator):
    _one_ion = True
    """
    Generate charge exchange spectra using the ACX2 model, for a single
    recombining ion. To use this model, you must have the acx2
    (https://acx2.readthedocs.io/) and pyatomdb (https://atomdb.readthedocs.io/)
    packages installed.

    Parameters
    ----------
    emin : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The minimum energy for the spectral model.
    emax : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
        The maximum energy for the spectral model.
    nbins : integer
        The number of bins in the spectral model.
    collntype : integer, optional
        The type of collision parameter to use:
        1 - center of mass energy (kev/u)
        2 - center of mass velocity (km/s)
        3 - donor ion velocity (km/s)
        4 - recombining ion velocity (km/s)
        Default: 1
    acx_model : integer, optional
        ACX model to fall back on, from 1 to 8. Default: 8
    recomb_type : integer, optional
        The type of recombination to use:
        1 - single recombination
        2 - full recombination
        3 - full renormalized recombination
        Default: 1
    binscale : string, optional
        The scale of the energy binning: "linear" or "log".
        Default: "linear"
    abund_table : string or array_like, optional
        The abundance table to be used for solar abundances.
        Either a string corresponding to a built-in table or an array
        of 30 floats corresponding to the abundances of each element
        relative to the abundance of H. Default is set in the SOXS
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
        "cl17.03" : the abundance table used in Cloudy v17.03.
    """

    def __init__(
        self,
        emin,
        emax,
        nbins,
        collntype=1,
        acx_model=8,
        recomb_type=1,
        binscale="linear",
        abund_table=None,
    ):
        super().__init__(
            emin,
            emax,
            nbins,
            collntype=collntype,
            acx_model=acx_model,
            recomb_type=recomb_type,
            binscale=binscale,
            var_elem=None,
            abund_table=abund_table,
        )

    def make_table(self, ions, colls, redshift):
        from astropy.units import Quantity

        numc = colls.size
        numi = len(ions)
        h_spec = np.zeros((numi, numc, self.nbins))
        he_spec = np.zeros((numi, numc, self.nbins))

        self.model.set_abund(np.ones(self.num_elements), elements=self.elements)

        # Shift the energy bins to the source frame
        self.model.set_ebins(self.ebins * (1.0 + redshift))

        # collision parameter parsing
        if not hasattr(colls, "unit"):
            colls = Quantity(colls, self.coll_units)
        colls = colls.to_value(self.coll_units)

        for i, (Z, ion) in enumerate(ions):
            mylog.info("Creating spectrum table for %s %i", elem_names[Z], ion)
            # Set the ionization fraction
            ionfrac = {}
            for e in self.model.elements:
                ionfrac[e] = np.zeros(e + 1)
            ionfrac[Z][ion] = 1.0
            self.model.set_ionfrac(ionfrac)
            pbar = tqdm(leave=True, total=numc, desc="Preparing spectrum table ")
            for j, coll in enumerate(colls):
                self.model.set_donorabund(["H", "He"], [1.0, 0.0])
                h_spec[i, j, :] = self.model.calc_spectrum(coll)
                self.model.set_donorabund(["H", "He"], [0.0, 1.0])
                he_spec[i, j, :] = self.model.calc_spectrum(coll)
                pbar.update()
            pbar.close()
        return h_spec, he_spec

    def get_spectrum(
        self, elem, ion, collnpar, He_frac, redshift, norm, velocity=0.0, tbroad=0.0
    ):
        """
        Get a charge exchange spectrum for a single recombining ion.

        Parameters
        ----------
        elem : string or integer
            The number or name of the element.
        ion : integer
            The ionization state of the element.
        collnpar : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`
            The collision parameter. Units are determined by the value of
            the collntype parameter set in the class constructor.
        He_frac : float
            Number fraction of donor which is He (remainder is H).
        redshift : float
            The redshift.
        norm : float
            The normalization of the model, in units of
            EM/(4*pi*(1+z)**2*D_A**2), where EM = int N_H^r N_(H+He)^d dV.
        velocity : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The velocity broadening parameter, in units of
            km/s. Default: 0.0
        tbroad : float, (value, unit) tuple, or :class:`~astropy.units.Quantity`, optional
            The temperature broadening parameter, in units of keV. Default: 0.0
        """
        # Get the atomic number
        if isinstance(elem, str):
            Z = elem_names.index(elem)
        else:
            Z = elem

        # check velocity and tbroad
        velocity = parse_value(velocity, "km/s")
        tbroad = parse_value(tbroad, "keV")

        # Set the ionization fraction
        ionfrac = {}
        for i in self.model.elements:
            ionfrac[i] = np.zeros(i + 1)
        ionfrac[Z][ion] = 1.0
        self.model.set_ionfrac(ionfrac)

        self.model.set_abund(np.ones(self.num_elements), elements=self.elements)

        spec = self._get_spectrum(redshift, He_frac, collnpar, tbroad, velocity)

        return Spectrum(self.ebins, norm * spec / self.de, binscale=self.binscale)
