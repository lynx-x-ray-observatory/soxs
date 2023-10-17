import soxs

spec = {
    "bi": soxs.Spectrum.from_xspec_script("aciss_stowed_model.xcm", 0.02, 15.0, 10000),
    "fi": soxs.Spectrum.from_xspec_script("acis-i-part-bg.xcm", 0.02, 15.0, 10000),
}
exp_time = (1.0, "Ms")

for cy in [0, 22]:
    rmf = "aciss_aimpt_cy%d.rmf" % cy
    instrument = (
        None,
        rmf,
    )
    for key in spec:
        out_file = rmf.replace(".rmf", f"_{key}.pha")
        soxs.simulate_spectrum(
            spec[key], instrument, exp_time, out_file, overwrite=True, noisy=False
        )
