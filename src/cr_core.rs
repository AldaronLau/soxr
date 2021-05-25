#[repr(C)]
struct HalfFirInfo {
    num_coefs: i32,
    coefs: *const f32,
    att: f32,
}

#[repr(C)]
struct CrCore {
    half_firs: *const HalfFirInfo,
    half_firs_len: usize,
    poly_firs: *const PolyFir,
}

#[repr(C)]
struct PolyFir {
    beta: f32,
    interp: PolyFir1,
}

#[repr(C)]
struct PolyFir1 {
    scalar: f32,
    stage_fn_t fn_: fn(input: *mut Stage, output: *mut Fifo),
}

#[repr(C)]
struct RateShared { /* So generated filter coefs may be shared between channels */
    poly_fir_coefs: *mut f32,
    dft_filter_t dft_filter[2];
}

#[repr(C)]
struct Stage {
    num: i32,

    /* Common to all stage types: */
    core_flags: i32,
    fn_: fn(input: *mut Stage, output: *mut Fifo),
    fifo: Fifo,
    pre: i32,       /* Number of past samples to store */
    pre_post: i32,  /* pre + number of future samples to store */
    preload: i32,   /* Number of zero samples to pre-load the fifo */
    out_in_ratio: f64, /* For buffer management. */
    input_size: i32,
    is_input: bool,

    /* For a stage with variable (run-time generated) filter coefs: */
    shared: *mut RateShared,
    dft_filter_num: u32, /* Which, if any, of the 2 DFT filters to use */
    dft_scratch: *mut f32,
    dft_out: *mut f32,
    coefs: *const f32,

    /* For a stage with variable L/M: */
    at: Step,
    step: Step,
    L: i32,
    remM: i32,
    n: i32,
    phase_bits: i32,
    block_len: i32,
    mult: f64,
    phase0: f64,
}

#[repr(C)]
struct Step { /* Fixed point arithmetic */
    ls: u64,
    ms: i64,
}

#[repr(C)]
struct Fifo {
    data: *mut u8,
    allocation: usize,   /* Number of bytes allocated for data. */
    item_size: usize,    /* Size of each item in data */
    begin: usize,        /* Offset of the first byte to read. */
    end: usize,          /* 1 + Offset of the last byte byte to read. */
}
