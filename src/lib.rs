use pyo3::prelude::*;
use rayon::prelude::*;

#[inline]
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            other => other,  // handle N, etc.
        })
        .collect()
}

/// Calculate the longest stop-free run in a single reading frame
#[inline]
fn longest_stop_free_in_frame(seq: &str, frame: usize) -> usize {
    let stop_codons = ["TAA", "TAG", "TGA"];
    let mut max_length = 0;
    let mut region_start = frame;
    let mut i = frame;

    while i + 3 <= seq.len() {
        let codon = &seq[i..i + 3];

        if stop_codons.contains(&codon) {
            let region_length = i - region_start;
            max_length = max_length.max(region_length);
            region_start = i + 3;
        }

        i += 3;
    }

    // Final region (after last stop, or entire frame if no stops)
    let final_region_length = i - region_start;
    max_length.max(final_region_length)
}

#[inline]
fn longest_stop_free_in_sequence(seq: &str) -> usize {
    (0..3)
        .map(|frame| longest_stop_free_in_frame(seq, frame))
        .max()
        .unwrap_or(0)
}

fn calculate_stop_free_run_single(seq: &str) ->usize {
    let seq_upper: String = seq.to_uppercase();
    let rc = reverse_complement(&seq_upper);
    
    let fwd_max = longest_stop_free_in_sequence(&seq_upper);
    let rc_max = longest_stop_free_in_sequence(&rc);

    fwd_max.max(rc_max)
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



