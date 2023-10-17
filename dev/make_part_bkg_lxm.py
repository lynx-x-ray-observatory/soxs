import soxs

spec = soxs.Spectrum.from_constant(3.0e-4, 0.02, 15.0, 10000)

exp_time = (1.0, "Ms")

for rmf in ["xrs_mucal_3.0eV.rmf", "xrs_mucal_1.5eV.rmf", "xrs_mucal_0.3eV.rmf"]:
    instrument = (
        None,
        rmf,
    )
    out_file = rmf.replace(".rmf", ".pha")

    soxs.simulate_spectrum(
        spec, instrument, exp_time, out_file, overwrite=True, noisy=False
    )
