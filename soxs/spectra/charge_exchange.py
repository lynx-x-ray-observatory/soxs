import numpy as np

from soxs.constants import abund_tables, cosmic_elem, elem_names, metal_elem
from soxs.utils import parse_value, soxs_cfg


class ChargeExchangeGenerator:
    def __init__(
        self,
        model,
        emin,
        emax,
        nbins,
        binscale="linear",
        var_elem=None,
        abund_table=None,
    ):
        try:
            from acx2 import ACXModel
        except ImportError:
            raise ImportError(
                "You must have the acx2 and pyatomdb packages "
                "installed to use the ChargeExchangeGenerator class!"
            )
        if model == "acx2":
            self.model = ACXModel()
        else:
            raise NotImplementedError(f"Model {model} not recognized!")
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
            self.var_elem = np.empty((0, 1), dtype="int")
        else:
            if len(var_elem) != len(set(var_elem)):
                raise RuntimeError(
                    'Duplicates were found in the "var_elem" list! %s' % var_elem
                )
            self.var_elem = [elem_names.index(e) for e in var_elem]
            self.var_elem.sort()
            self.var_elem = np.array(self.var_elem, dtype="int")
            self.var_elem_names = [elem_names[e[0]] for e in self.var_elem]
        self.num_var_elem = len(self.var_elem)
        self.cosmic_elem = [elem for elem in cosmic_elem if elem not in self.var_elem]
        self.metal_elem = [elem for elem in metal_elem if elem not in self.var_elem]
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

    def get_spectrum(self, kT, abund):
        pass
