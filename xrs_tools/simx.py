from __future__ import print_function
import os

def run_simx(simput_file, output_prefix, exposure_time, sky_center,
             instrument='hdxi', bkgnd_scale=1.0, random_seed=None):
    """
    Run SIMX assuming X-ray Surveyor mission parameters.

    Parameters
    ----------
    simput_file : string
        The path to the SIMPUT file used for the input. 
    output_prefix : string
        The prefix of the events file which will be written. The
        events file will be called *output_prefix*_evt.fits.
    exposure_time : float
        The exposure time in seconds. 
    sky_center : array-like
        The RA and Dec of the center pointing in degrees. Can be
        a tuple, list, or NumPy array.
    instrument : string, optional
        The X-ray Surveyor detector to use. Can be one of: "hdxi",
        "xcal", "cat". Case-insensitive. Default: "hdxi".
    bkgnd_scale : float, optional
        The scale factor to multiply the instrumental background by.
        Default: 1.0
    random_seed : integer, optional
        Set this to specify a random seed to use when drawing the 
        photons. Default None, which means that the random seed
        will depend on the system time. 
    """
    if instrument.lower() not in ["hdxi", "xcal", "cat"]:
        raise RuntimeError("Instrument %s not valid!" % instrument)
    os.system('punlearn simx')
    os.system('pset simx mode=hl')
    os.system("pset simx Exposure=%g" % exposure_time)
    os.system("pset simx UseSimput=yes")
    os.system("pset simx MissionName=XraySurveyor InstrumentName=%s" % instrument.upper())
    os.system("pset simx ScaleBkgnd=%g" % bkgnd_scale)
    if random_seed is not None:
        os.system("pset simx RandomSeed=%d" % random_seed)
    os.system("pset simx SimputFile=%s" % simput_file)
    os.system("pset simx PointingRA=%g" % sky_center[0])
    os.system("pset simx PointingDec=%g" % sky_center[1])
    os.system("pset simx OutputFileName=%s" % output_prefix)
    os.system("simx")
