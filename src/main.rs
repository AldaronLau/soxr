use std::os::raw::{
    c_void, c_char, c_uint
};

extern "C" {
    /* Create a stream resampler: */
    fn soxr_create(
        input_rate: f64,      /* Input sample-rate. */
        output_rate: f64,     /* Output sample-rate. */
        num_channels: c_uint, /* Number of channels to be used. */
        /* All following arguments are optional (may be set to NULL). */
        error: *mut *const c_char,   /* To report any error during creation. */
        spec: *const c_void,      /* To specify non-default I/O formats. */
        quality_spec: *const c_void, /* To specify non-default resampling quality.*/
        runtime_spec: *const c_void,  /* To specify non-default runtime resources. */
    ) -> *mut c_void;               
/*
        Default io_spec      is per soxr_io_spec(SOXR_FLOAT32_I, SOXR_FLOAT32_I)
        Default quality_spec is per soxr_quality_spec(SOXR_HQ, 0)
        Default runtime_spec is per soxr_runtime_spec(1)                          */

    /* If not using an app-supplied input function, after creating a stream
     * resampler, repeatedly call: */

    fn soxr_process(
        resampler: *mut c_void,      /* As returned by soxr_create. */
                                /* Input (to be resampled): */
        in_: *const c_void,             /* Input buffer(s); may be NULL (see below). */
        ilen: usize,           /* Input buf. length (samples per channel). */
        idone: *mut usize,        /* To return actual # samples used (<= ilen). */
                                /* Output (resampled): */
        out: *mut c_void,            /* Output buffer(s).*/
        olen: usize,           /* Output buf. length (samples per channel). */
        odone: *mut usize) -> *const c_char;       /* To return actual # samples out (<= olen).

        Note that no special meaning is associated with ilen or olen equal to
        zero.  End-of-input (i.e. no data is available nor shall be available)
        may be indicated by seting `in' to NULL.                                  */
}

fn resample_channel(input: Vec<f32>) -> Vec<f32> {
    let out_size = input.len() as f64 * (44_100.0 / 48_000.0);
    let mut output = vec![0.0; out_size as usize];

    // Sample to 44100 from 48000
    let resampler = unsafe {
        soxr_create(
            48_000.0,
            44_100.0,
            1, // Always going to be one.
            std::ptr::null_mut(),
            std::ptr::null(),
            std::ptr::null(),
            std::ptr::null(),
        )
    };

    let mut ilen = 0;
    let mut olen = 0;

    // Resample the whole thing at once.
    unsafe {
        // Add input
        soxr_process(
            resampler,
            input.as_ptr().cast(),
            input.len(),
            &mut ilen,
            output.as_mut_ptr().cast(),
            output.len(),
            &mut olen,
        );
        let mut flush_len = 0;
        // Flush output
        soxr_process(
            resampler,
            std::ptr::null(),
            0,
            &mut 0,
            output.as_mut_ptr().cast(),
            output.len(),
            &mut flush_len,
        );
        olen += flush_len;
    }

    // Should have used all of both input and output
    assert_eq!(ilen, input.len());
    assert_eq!(olen, output.len());
    
    output
}

fn main() {
    // Load 48_000 HZ song S16 Interleaved Stereo.
    let song = std::fs::read("AUSTRA - Reconcile.raw").expect("No file!");

    // First, de-interleave and convert to f32 (left and right channels).
    println!("Preparing…");

    let mut left = Vec::new();
    let mut right = Vec::new();

    for frame in song.chunks(4) {
        let l = i16::from_le_bytes([frame[0], frame[1]]);
        let r = i16::from_le_bytes([frame[2], frame[3]]);
        
        let l = (l as f32 + 0.5) / (i16::MAX as f32 + 0.5);
        let r = (r as f32 + 0.5) / (i16::MAX as f32 + 0.5);
        
        left.push(l);
        right.push(r);
    }

    // Process on separate threads.
    println!("Resampling…");
    
    let left = std::thread::spawn(|| resample_channel(left));
    let right = std::thread::spawn(|| resample_channel(right));

    // Get the resampled left and right channels individually.
    let left: Vec<f32> = left.join().unwrap();
    let right: Vec<f32> = right.join().unwrap();

    println!("Resampling Finished!");

    let mut out: Vec<u8> = Vec::new();

    // Interleave samples.
    for (l, r) in left.iter().zip(right.iter()) {
        let l = ((l * (i16::MAX as f32 + 0.5)) - 0.5) as i16;
        let r = ((r * (i16::MAX as f32 + 0.5)) - 0.5) as i16;

        out.extend(&l.to_le_bytes());
        out.extend(&r.to_le_bytes());
    }

    // Reading test....
    let song = std::fs::read("test.raw").expect("No test!");

    assert_eq!(song, out);

/*
    println!("Writing out...");

    std::fs::write("out.raw", out).expect("Could not write!");*/
}
