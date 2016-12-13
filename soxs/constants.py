import astropy.units as apu
from astropy.constants import h, c, u
import numpy as np

one_arcsec = 1.0/3600.0

erg_per_eV = apu.eV.to("erg")
erg_per_keV = erg_per_eV * 1.0e3
keV_per_erg = 1.0 / erg_per_keV
eV_per_erg = 1.0 / erg_per_eV

hc = (h*c).to("keV*angstrom").value
clight = c.to("cm/s").value

m_u = u.to("g").value

cosmic_elem = [1,2,3,4,5,9,11,15,17,19,21,22,23,24,25,27,29,30]
metal_elem = [6,7,8,10,12,13,14,16,18,20,26,28]

atomic_weights = np.array([0.0,1.00794,4.00262,6.941,9.012182,10.811,
                           12.0107,14.0067,15.9994,18.9984,20.1797,
                           22.9898,24.3050,26.9815,28.0855,30.9738,
                           32.0650,35.4530,39.9480,39.0983,40.0780,
                           44.9559,47.8670,50.9415,51.9961,54.9380,
                           55.8450,58.9332,58.6934,63.5460,65.3800])
