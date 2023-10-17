import soxs

spec = soxs.Spectrum.from_xspec_script("acis-i-part-bg.xcm", 0.02, 15.0, 10000)

exp_time = (1.0, "Ms")

for cy in [0, 22]:
    rmf = "acisi_aimpt_cy%d.rmf" % cy
    instrument = (
        None,
        rmf,
    )
    out_file = rmf.replace(".rmf", ".pha")

    soxs.simulate_spectrum(
        spec, instrument, exp_time, out_file, overwrite=True, noisy=False
    )
