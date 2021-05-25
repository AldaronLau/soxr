fn main() {
    cc::Build::new()
        .file("soxr/cr.c")
        .file("soxr/cr32.c")
        .file("soxr/soxr.c")
        .file("soxr/dbesi0.c")
        .file("soxr/filter.c")
        .file("soxr/fft4g32.c")
        .warnings(true)
        .extra_warnings(true)
        .flag("-Wall")
/*
        .warnings(false)
        .extra_warnings(false)
        .flag("-w")
*/
        .compile("soxr");
}
