from __future__ import print_function
import os

def run_simx(simput_file, output_prefix, exposure_time, sky_center,
             instrument='hdxi', bkgnd_scale=1.0, random_seed=None):
    if instrument.lower() not in ["hdxi", "xcal", "cat"]:
        raise RuntimeError("Instrument %s not valid!" % instrument)
    os.system('punlearn simx')
    os.system('pset simx mode=hl')
    os.system("pset simx Exposure=%g" % exposure_time)
    os.system("pset simx UseSimput=yes")
    os.system("pset simx MissionName=XraySurveyor InstrumentName=%g" % instrument.upper())
    os.system("pset simx ScaleBkgnd=%g" % bkgnd_scale)
    if random_seed is not None:
        os.system("pset simx RandomSeed=%d" % random_seed)
    os.system("pset simx SimputFile=%s" % simput_file)
    os.system("pset simx PointingRA=%g PointingDec=%g" % sky_center)
    os.system("pset simx OutputFileName=%s" % output_prefix)
    os.system("simx")
