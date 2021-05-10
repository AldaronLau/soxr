fn main() {
    cc::Build::new()
        .file("soxr/cr.c")
        .file("soxr/cr32.c")
        .file("soxr/soxr.c")
        .file("soxr/dbesi0.c")
        .file("soxr/filter.c")
        .file("soxr/data-io.c")
        .file("soxr/fft4g32.c")
        .file("soxr/fft4g64.c")
        .file("soxr/pffft-wrap.c")
        // These are included by other C Files.
        // SIMD!
//.file("soxr/fft4g.c")
//.file("soxr/pffft.c")
//.file("soxr/cr-core.c")
//.file("soxr/ui")
        .warnings(false)
        .extra_warnings(false)
        .flag("-w")
        .define("HAVE_BIGENDIAN", Some("0"))
        .compile("soxr");
}
