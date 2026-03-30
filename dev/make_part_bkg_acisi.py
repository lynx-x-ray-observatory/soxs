import argparse

import soxs

parser = argparse.ArgumentParser()
parser.add_argument("cycles", type=str)
args = parser.parse_args()

cycles = [int(c.strip()) for c in args.cycles.split(",")]

spec = soxs.Spectrum.from_xspec_script("acis-i-part-bg.xcm", 0.02, 15.0, 10000)

exp_time = (1.0, "Ms")

for cy in cycles:
    rmf = f"acisi_aimpt_cy{cy}.rmf"
    instrument = (
        None,
        rmf,
    )
    out_file = rmf.replace(".rmf", ".pha")

    soxs.simulate_spectrum(spec, instrument, exp_time, out_file, overwrite=True, noisy=False)
