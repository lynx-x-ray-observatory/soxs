import numpy as np
from tqdm.auto import tqdm

from soxs.constants import abund_tables, elem_names
from soxs.spectra import Spectrum
from soxs.utils import parse_value, soxs_cfg

elem_subset = [1, 2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28]
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
        if var_elem is None:
            self.var_elem = np.empty(0, dtype="int")
            self.elements = elem_full if self._one_ion else elem_subset
        else:
            if len(var_elem) != len(set(var_elem)):
                raise RuntimeError(
                    'Duplicates were found in the "var_elem" list! %s' % var_elem
                )
            self.var_elem = [elem_names.index(e) for e in var_elem]
            self.var_elem.sort()
            self.var_elem = np.array(self.var_elem, dtype="int")
            self.var_elem_names = [elem_names[e] for e in self.var_elem]
            self.elements = elem_full
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
            self._cv = np.sqrt(4786031.3 / 25.0)
        else:
            self.coll_units = "km/s"
            if collntype == 2:
                self._cv = 1.0
            elif collntype == 3:
                self._cv = 1.0 / 13.0
            elif collntype == 4:
                self._cv = 12.0 / 13.0

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

        # normalize the spectrum
        cv = collnpar
        if self.collntype == 1:
            cv **= 0.5
        cv *= self._cv

        return spec / cv

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
        numc = colls.size
        numi = ions.size
        h_spec = np.zeros((numi, numc, self.nbins))
        he_spec = np.zeros((numi, numc, self.nbins))

        self.model.set_abund(np.ones(self.num_elements), elements=self.elements)

        pbar = tqdm(leave=True, total=numc, desc="Preparing spectrum table ")
        for i, (Z, ion) in enumerate(ions):
            # Set the ionization fraction
            ionfrac = {}
            for i in self.model.elements:
                ionfrac[i] = np.zeros(i + 1)
            ionfrac[Z][ion] = 1.0
            self.model.set_ionfrac(ionfrac)
            for j, collnpar in enumerate(colls):
                h_spec[i, j, :] = self._get_spectrum(redshift, 0.0, collnpar, 0.0, 0.0)
                he_spec[i, j, :] = self._get_spectrum(redshift, 1.0, collnpar, 0.0, 0.0)
            pbar.update()
        pbar.close()
        return h_spec, he_spec

    def get_spectrum(
        self, elem, ion, collnpar, He_frac, redshift, norm, velocity=0.0, tbroad=0.0
    ):

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
