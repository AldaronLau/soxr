use std::os::raw::{c_void};

extern "C" {
    fn resampler_input(rate: *mut c_void, samples: *const f32, n: usize) -> *mut f32;
    fn resampler_process(rate: *mut c_void, olen: usize);
    fn resampler_output(rate: *mut c_void, samples: *mut f32, n0: *mut usize) -> *const f32;
    fn resampler_flush(rate: *mut c_void);
}

/// Process audio adding input and writing to output as many samples as possible.
pub(crate) fn process(resampler: *mut c_void, in_: &[f32], out: &mut [f32]) -> usize {
    let idone = input(resampler, in_.as_ptr(), in_.len());
    let odone = output(resampler, out.as_mut_ptr(), out.len());

    assert_eq!(idone, in_.len());

    odone
}

/// Flush to output.
pub(crate) fn flush(resampler: *mut c_void, out: &mut [f32]) {
    unsafe {
        resampler_flush(resampler);
    }

    let odone = output(resampler, out.as_mut_ptr(), out.len());
    
    assert_eq!(odone, out.len());
}

fn input(resampler: *mut c_void, inp: *const f32, len: usize) -> usize {
    unsafe {
        let chan = resampler_input(resampler, std::ptr::null(), len);
        std::ptr::copy_nonoverlapping(inp, chan, len);
    };

    return len;
}

fn output(resampler: *mut c_void, out: *mut f32, mut len: usize) -> usize {
    unsafe {
        resampler_process(resampler, len);
        let src = resampler_output(resampler, std::ptr::null_mut(), &mut len);
        std::ptr::copy_nonoverlapping(src, out, len);
    }

    return len;
}
