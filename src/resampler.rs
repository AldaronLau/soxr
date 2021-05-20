use std::os::raw::{c_void};

extern "C" {
    fn resampler_input(rate: *mut c_void, samples: *const f32, n: usize) -> *mut f32;
    fn resampler_process(rate: *mut c_void, olen: usize);
    fn resampler_output(rate: *mut c_void, samples: *mut f32, n0: *mut usize) -> *const f32;
    fn resampler_flush(rate: *mut c_void);
}

/// Process audio adding input and writing to output as many samples as possible.
pub(crate) fn process(resampler: *mut c_void, in_: &[f32], out: &mut [f32]) -> usize {
    input(resampler, in_);
    output(resampler, out)
}

/// Flush to output.
pub(crate) fn flush(resampler: *mut c_void, out: &mut [f32]) -> usize {
    unsafe {
        resampler_flush(resampler);
    }

    output(resampler, out)
}

fn input(resampler: *mut c_void, inp: &[f32]) {
    unsafe {
        resampler_input(resampler, inp.as_ptr(), inp.len());
    };
}

fn output(resampler: *mut c_void, out: &mut [f32]) -> usize {
    let mut len = out.len();

    unsafe {
        resampler_process(resampler, out.len());
        resampler_output(resampler, out.as_mut_ptr(), &mut len);
    }

    return len;
}
