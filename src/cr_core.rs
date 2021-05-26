use std::os::raw::{c_void, c_char};

#[repr(C)]
struct HalfFirInfo {
    num_coefs: i32,
    coefs: *const f32,
    att: f32,
}

#[repr(C)]
pub(crate) struct CrCore {
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
    fn_: unsafe extern "C" fn(input: *mut Stage, output: *mut Fifo),
}

#[repr(C)]
struct DftFilter {
    dft_length: i32,
    num_taps: i32,
    post_peak: i32,
    coefs: *mut f32,
}

#[repr(C)]
pub(crate) struct RateShared { /* So generated filter coefs may be shared between channels */
    poly_fir_coefs: *mut f32,
    dft_filter: [DftFilter; 2],
}

#[repr(C)]
pub(crate) struct Stage {
    num: i32,

    /* Common to all stage types: */
    core_flags: i32,
    fn_: unsafe extern "C" fn(input: *mut Stage, output: *mut Fifo),
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
    l: i32,
    rem_m: i32,
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

//// Constants ////

const HALF_FIR_COEFS_7: &[f32] = &[
     3.1062656496657370e-01, -8.4998810699955796e-02,  3.4007044621123500e-02,
    -1.2839903789829387e-02,  3.9899380181723145e-03, -8.9355202017945374e-04,
     1.0918292424806546e-04,
];

const HALF_FIR_COEFS_8: &[f32] = &[
     3.1154652365332069e-01, -8.7344917685739543e-02,  3.6814458353637280e-02,
    -1.5189204581464479e-02,  5.4540855610738801e-03, -1.5643862626630416e-03,
     3.1816575906323303e-04, -3.4799449225005688e-05,
];

const HALF_FIR_COEFS_9: &[f32] = &[
     3.1227034755311189e-01, -8.9221517147969526e-02,  3.9139704015071934e-02,
    -1.7250558515852023e-02,  6.8589440230476112e-03, -2.3045049636430419e-03,
     6.0963740543348963e-04, -1.1323803957431231e-04,  1.1197769991000046e-05,
];

const HALF_FIRS: &[HalfFirInfo] = &[
    HalfFirInfo {
        num_coefs: 7,
        coefs: HALF_FIR_COEFS_7.as_ptr(),
        att: 120.65,
    },
    HalfFirInfo {
        num_coefs: 8,
        coefs: HALF_FIR_COEFS_8.as_ptr(),
        att: 136.51,
    },
    HalfFirInfo {
        num_coefs: 9,
        coefs: HALF_FIR_COEFS_9.as_ptr(),
        att: 152.32,
    }
];

unsafe fn fifo_occupancy(f: *mut Fifo) -> i32 {
    (((*f).end - (*f).begin) / (*f).item_size) as i32
}

unsafe fn fifo_clear(f: *mut Fifo) {
    (*f).begin = 0;
    (*f).end = 0;
}

unsafe fn fifo_reserve(f: *mut Fifo, n0: i32) -> *mut c_void {
    let n = n0 as usize * (*f).item_size;

    if (*f).begin == (*f).end {
        fifo_clear(f);
    }

    loop {
        if (*f).end + n <= (*f).allocation {
            let p: *mut c_void = (*f).data.offset((*f).end as isize).cast();

            (*f).end += n;
            return p;
        }
        if (*f).begin > 0x4000 /* fifo min */ {
            memmove(
                (*f).data.cast(),
                (*f).data.offset((*f).begin as isize).cast(),
                (*f).end - (*f).begin
            );
            (*f).end -= (*f).begin;
            (*f).begin = 0;
            continue;
        }
        (*f).data = realloc((*f).data.cast(), (*f).allocation + n).cast();
        (*f).allocation += n;
        if (*f).data.is_null() {
            return std::ptr::null_mut();
        }
    }
}

unsafe fn fifo_read(f: *mut Fifo, n0: i32, data: *mut c_void) -> *mut c_void {
    let n = n0 as usize * (*f).item_size;
    let ret: *mut c_char = (*f).data.offset((*f).begin as isize).cast();
    if n > ((*f).end - (*f).begin) {
        return std::ptr::null_mut();
    }
    if !data.is_null() {
        memcpy(data, ret.cast(), n as usize);
    }
    (*f).begin += n;
    return ret.cast();
}

unsafe extern "C" fn vpoly0(p: *mut Stage, output_fifo: *mut Fifo) {
    let fifo_occupancy = fifo_occupancy(&mut (*p).fifo);
    let stage_occupancy = (fifo_occupancy - (*p).pre_post).max(0);
    let num_in: i32 = stage_occupancy.min((*p).input_size);
    
    if num_in != 0 {
        let input: *const f32 = fifo_read(&mut (*p).fifo, 0, std::ptr::null_mut()).offset((*p).pre as isize).cast();
        let mut at: i32 = ((*p).at.ms >> 32) as i32;
        let step: i32 = ((*p).step.ms >> 32) as i32;
        let num_out: i32 = (num_in * (*p).l - at + step - 1) / step;
        let output: *mut f32 = fifo_reserve(output_fifo, num_out).cast();

        let mut i = 0;
        while at < num_in * (*p).l {
            let div = at / (*p).l;
            let rem = at % (*p).l;
            let at2: *const f32 = input.offset(div as isize);
            let mut j = 0;
            let mut sum = 0.0;
            let coefs: *mut f32 = (*(*p).shared).poly_fir_coefs.offset(((*p).n * rem) as isize);
            while j < ((*p).n) {
                sum += *coefs.offset(j as isize) * *at2.offset(j as isize);
                j += 1;
            }
            *output.offset(i as isize) = sum;
            i += 1;
            at += step;
        }

        assert_eq!(i, num_out);
        fifo_read(&mut (*p).fifo, at / (*p).l, std::ptr::null_mut());
        (*p).at.ms = ((*p).at.ms & 0x00000000FFFFFFFF)
            | ((at % (*p).l) as i64) << 32;
    }
}

unsafe extern "C" fn u100_0(p: *mut Stage, output_fifo: *mut Fifo) {
    let fifo_occupancy = fifo_occupancy(&mut (*p).fifo);
    let stage_occupancy = (fifo_occupancy - (*p).pre_post).max(0);
    let num_in: i32 = stage_occupancy.min((*p).input_size);
    
    if num_in != 0 {
        let mut at = ((*p).at.ms >> 32) as i32;
        let step = ((*p).step.ms >> 32) as i32;
        let num_out = ((num_in * (*p).l - at + step - 1) / step) as i32;
        let output: *mut f32 = fifo_reserve(output_fifo, num_out).cast();

        let mut i = 0;
        while at < num_in * (*p).l {
            *output.offset(i) = 0.0;
            i += 1;
            at += step;
        }

        assert_eq!(i, num_out as isize);
        fifo_read(&mut (*p).fifo, at / (*p).l, std::ptr::null_mut());
        (*p).at.ms = ((*p).at.ms & 0x00000000FFFFFFFF)
            | ((at % (*p).l) as i64) << 32;
    }
}

const POLY_FIRS: &[PolyFir] = &[
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0 } },
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},

    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},

    PolyFir { beta: 10.62, interp: PolyFir1 { scalar: 42.0, fn_: u100_0}},
    PolyFir { beta: 11.28, interp: PolyFir1 { scalar: 11.0, fn_: u100_0}},

    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
    PolyFir { beta: -1.0, interp: PolyFir1 { scalar: 0.0, fn_: vpoly0}},
];

// 32-bit Floating Point
pub(crate) const CR32: &CrCore = &CrCore {
    half_firs: HALF_FIRS.as_ptr(),
    half_firs_len: HALF_FIRS.len(),
    poly_firs: POLY_FIRS.as_ptr(),
};

//// Functions

extern "C" {
    fn realloc(ptr: *mut c_void, size: usize) -> *mut c_void;
    fn memmove(dest: *mut c_void, src: *const c_void, n: usize) -> *mut c_void;
    fn memcpy(dest: *mut c_void, src: *const c_void, n: usize) -> *mut c_void;
}
