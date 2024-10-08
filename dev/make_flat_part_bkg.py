import argparse

import soxs

parser = argparse.ArgumentParser(description="Make a flat particle background spectrum")
parser.add_argument(
    "const_flux",
    type=float,
    help="The constant flux in photons/s/keV, assuming 1 arcmin^2 solid angle.",
)
parser.add_argument("rmf", type=str, help="The RMF to use.")
parser.add_argument("outfile", type=str, help="The output file name.")
args = parser.parse_args()

spec = soxs.Spectrum.from_constant(args.const_flux, 0.01, 20.0, 60000)

exp_time = (1.0, "Ms")

instrument = (
    None,
    args.rmf,
)

soxs.simulate_spectrum(
    spec, instrument, exp_time, args.outfile, overwrite=True, noisy=False
)
