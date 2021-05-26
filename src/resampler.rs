use std::os::raw::{c_char};

use crate::cr_core::{CR32, CrCore, Stage, RateShared};

#[repr(C)]
pub(super) struct Resampler {
    core: *const CrCore,
    io_ratio: f64,
    samples_in: i64,
    samples_out: i64,
    num_stages: i32,
    flushing: i32,
    stages: *mut Stage,

    shared: RateShared, /* Between channels. */
}

extern "C" {
    fn resampler_input(rate: *mut Resampler, samples: *const f32, n: usize)
        -> *mut f32;
    fn resampler_process(rate: *mut Resampler, olen: usize);
    fn resampler_output(rate: *mut Resampler, samples: *mut f32, n0: *mut usize)
        -> *const f32;
    fn resampler_flush(rate: *mut Resampler);
    fn resampler_init(
        p: *mut Resampler,         /* Per audio channel. */
        shared: *mut RateShared,    /* By channels undergoing same rate change. */
        io_ratio: f64,          /* Input rate divided by output rate. */
        core: *const CrCore,
    ) -> *const c_char;
}

pub(crate) fn new(io_ratio: f64) -> *mut Resampler {
    unsafe {
        let resampler = Box::<Resampler>::new(std::mem::zeroed());
        let resampler = Box::leak(resampler);

        resampler_init(resampler, &mut (*resampler).shared, io_ratio, CR32);

        resampler
    }
}

/// Process audio adding input and writing to output as many samples as possible.
pub(crate) fn process(resampler: *mut Resampler, in_: &[f32], out: &mut [f32]) -> usize {
    input(resampler, in_);
    output(resampler, out)
}

/// Flush to output.
pub(crate) fn flush(resampler: *mut Resampler, out: &mut [f32]) -> usize {
    unsafe {
        resampler_flush(resampler);
    }

    output(resampler, out)
}

fn input(resampler: *mut Resampler, inp: &[f32]) {
    unsafe {
        resampler_input(resampler, inp.as_ptr(), inp.len());
    };
}

fn output(resampler: *mut Resampler, out: &mut [f32]) -> usize {
    let mut len = out.len();

    unsafe {
        resampler_process(resampler, out.len());
        resampler_output(resampler, out.as_mut_ptr(), &mut len);
    }

    return len;
}
