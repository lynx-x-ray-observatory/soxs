Spectrum API
============

.. autoclass:: soxs.spectra.base.Spectrum
    :members:

.. autoclass:: soxs.spectra.base.CountRateSpectrum
    :members: generate_energies, rescale_flux, new_spec_from_band, get_lum_in_band, from_spectrum, from_powerlaw, from_constant

.. autoclass:: soxs.spectra.base.ConvolvedSpectrum
    :members: deconvolve, generate_energies, rescale_flux, new_spec_from_band
