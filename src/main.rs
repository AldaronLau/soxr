use std::thread;

mod resampler;
mod cr_core;

fn resample_channel(input: Vec<f32>, hz_in: f64, hz_out: f64) -> Vec<f32> {
    let out_size = input.len() as f64 * (hz_out / hz_in);
    let mut output = vec![0.0; out_size.round() as usize];

    // Sample to 44100 from 48000
    let mut resampler = resampler::new(hz_in / hz_out);

    // Resample the whole thing at once.
    let x = resampler::process(&mut resampler, &input, &mut output);
    resampler::flush(&mut resampler, &mut output[x..]);

    output
}

fn resample_audio(filename: &str, hz_in: f64, hz_out: f64) -> Vec<u8> {
    // Load song S16 Interleaved Stereo.
    let song = std::fs::read(filename).expect("No file!");

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
    
    let left = thread::spawn(move || resample_channel(left, hz_in, hz_out));
    let right = thread::spawn(move || resample_channel(right, hz_in, hz_out));

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

    out
}

fn main() {
    // Downsample
    let austra = resample_audio("AUSTRA - Reconcile.raw", 48_000.0, 44_100.0);
    // Upsample
    let shaed = resample_audio("SHAED - ISOU.raw", 44_100.0, 48_000.0);

    // Reading test....
    //std::fs::write("test.raw", &austra).expect("Could not write!");
    let song = std::fs::read("test.raw").expect("No test!");
    assert!(song == austra);

    // Reading check...
    //std::fs::write("check.raw", &shaed).expect("Could not write!");
    let song = std::fs::read("check.raw").expect("No test!");
    assert!(song == shaed);
}
