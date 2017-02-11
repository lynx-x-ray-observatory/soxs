from soxs.instrument import make_background, AuxiliaryResponseFile, \
    FlatResponse
from soxs.background.foreground import hm_astro_bkgnd
from soxs.background.instrument import acisi_particle_bkgnd
from soxs.background.spectra import ConvolvedBackgroundSpectrum
from numpy.random import RandomState
import numpy as np

prng = RandomState(29)

def test_uniform_bkgnd_scale():
    hdxi_arf = AuxiliaryResponseFile("xrs_hdxi_3x10.arf")
    flat_arf = FlatResponse(hdxi_arf.elo[0], hdxi_arf.ehi[-1], 1.0, 
                            hdxi_arf.emid.size)
    events, event_params = make_background(50000.0, "hdxi", [30., 45.], 
                                           foreground=True, instr_bkgnd=True,
                                           ptsrc_bkgnd=False, cosmo_bkgnd=False,
                                           prng=prng)
    ncts = np.logical_and(events["energy"] >= 0.7, events["energy"] <= 2.0).sum()
    t_exp = event_params["exposure_time"]
    fov = (event_params["fov"]*60.0)**2
    S = ncts/t_exp/fov
    dS = np.sqrt(ncts)/t_exp/fov
    foreground = ConvolvedBackgroundSpectrum(hm_astro_bkgnd, hdxi_arf)
    instr_bkgnd = ConvolvedBackgroundSpectrum(acisi_particle_bkgnd, flat_arf)
    f_sum = foreground.get_flux_in_band(0.7, 2.0)[0]
    i_sum = instr_bkgnd.get_flux_in_band(0.7, 2.0)[0]
    b_sum = (f_sum+i_sum).to("ph/(arcsec**2*s)").value
    assert np.abs(S-b_sum) < dS
    
if __name__ == "__main__":
    test_uniform_bkgnd_scale()