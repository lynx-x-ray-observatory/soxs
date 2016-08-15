import astropy.units as u

erg_per_eV = u.eV.to("erg")
erg_per_keV = erg_per_eV * 1.0e3
keV_per_erg = 1.0 / erg_per_keV
eV_per_erg = 1.0 / erg_per_eV
