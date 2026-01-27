use pyo3::prelude::*;
use rayon::prelude::*;
use regex::Regex;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref TAA_STOP: Regex = Regex::new(r"TAA").unwrap();
    pub static ref TAG_STOP: Regex = Regex::new(r"TAG").unwrap();
    pub static ref TGA_STOP: Regex = Regex::new(r"TGA").unwrap();
}

fn calculate_stop_free_run_single(seq: &str) ->usize {
    let rev_seq = seq.chars().rev().collect::<String>();

    let frames = [
        &seq[0..],
        &seq[1.min(seq.len())..],
        &seq[2.min(seq.len())..],
        &rev_seq[0..],
        &rev_seq[1.min(rev_seq.len())..],
        &rev_seq[2.min(rev_seq.len())..],
    ];

    let stops: Vec<&Regex> = vec![&TAA_STOP, &TAG_STOP, &TGA_STOP];
    frames
        .par_iter().map(|frame| 
            stops
            .iter()
            .map(|x| x.find(&frame))
            .filter_map(|x| x.map(|y| y.end()))
            .collect::<Vec<usize>>()
            .into_iter()
            .max()
            .unwrap_or(frame.len())
        ).collect::<Vec<usize>>()
        .into_iter()
        .max()
        .unwrap_or(seq.len())
}

#[pyfunction]
fn calculate_stop_free_run(sequences: Vec<String>) -> PyResult<Vec<usize>> {
    Ok(
        sequences
        .par_iter()
        .map(|seq| calculate_stop_free_run_single(seq))
        .collect()
    )
}

#[pyfunction]
fn calculate_stop_free_runs_with_ids(
    seq_tuples: Vec<(String, String)>
) -> PyResult<Vec<(String, usize)>> {
    Ok(seq_tuples
        .par_iter()
        .map(|(id, seq)| (id.clone(), calculate_stop_free_run_single(seq)))
        .collect())
}


/// A Rust implementation of stopFree
#[pymodule]
fn stopfree(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calculate_stop_free_run, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_stop_free_runs_with_ids, m)?)?;
    Ok(())
}



