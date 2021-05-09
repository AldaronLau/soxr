fn main() {
    cc::Build::new()
        .file("soxr/cr.c")
        .file("soxr/cr32.c")
        //.file("soxr/cr64.c")
        .file("soxr/soxr.c")
        .file("soxr/dbesi0.c")
        .file("soxr/filter.c")
        .file("soxr/data-io.c")
        .file("soxr/fft4g32.c")
        .file("soxr/fft4g64.c")
        .file("soxr/soxr-lsr.c")
        .file("soxr/pffft-wrap.c")
        // These are included by other C Files.
        // SIMD!
//.file("soxr/pffft32.c") // defines _soxr_rdft32_cb
//.file("soxr/fft4g.c")
//.file("soxr/pffft.c")
        //.file("soxr/avfft32.c")
        // This file requires WITH_VR32
//        .file("soxr/vr32.c")
        // This file requires WITH_CR32S
//        .file("soxr/cr32s.c")
        // This file requires WITH_CR64S
//        .file("soxr/cr64s.c")
        // This file requires AVCODEC_FOUND?

        // This file is included by other .c files
//        .file("soxr/cr-core.c")
        // This file requires SIMD
//        .file("soxr/util32s.c")
        // This file requires SIMD
//        .file("soxr/util64s.c")
        // This file requires SIMD?
//        .file("soxr/avfft32s.c")
        // This file requires SIMD
//        .file("soxr/fft4g32s.c")
        // This file requires SIMD
//        .file("soxr/pffft32s.c")
        // This file requires SIMD
//        .file("soxr/pffft64s.c")
        // This files is for generating variable rate coefs, not-needed
//        .file("soxr/vr-coefs.c")
        // This file requires SIMD
//        .file("soxr/util-simd.c")
        .warnings(false)
        .extra_warnings(false)
        .flag("-w")
        .define("AVCODEC_FOUND", Some("0"))
        .define("AVUTIL_FOUND", Some("0"))
        .define("WITH_PFFFT", Some("1"))
        .define("HAVE_FENV_H", Some("1"))
        .define("HAVE_LRINT", Some("1"))
        .define("HAVE_BIGENDIAN", Some("0"))
        .define("WITH_CR32", Some("1"))
        .define("WITH_CR32S", Some("0"))
        .define("WITH_CR64", Some("0"))
        .define("WITH_CR64S", Some("0"))
        .define("WITH_VR32", Some("0"))
        .define("WITH_HI_PREC_CLOCK", Some("1"))
        .define("WITH_FLOAT_STD_PREC_CLOCK", Some("0"))
        .define("WITH_DEV_TRACE", Some("0"))
        .define("SOXR_LIB", Some("1"))
        .compile("soxr");
}
