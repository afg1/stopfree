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

/// Calculate the longest stop-free run in a single reading frame in codons
#[inline]
fn longest_stop_free_in_frame(seq: &str, frame: usize) -> usize {
    let stop_codons = ["TAA", "TAG", "TGA"];
    let mut max_length = 0;
    let mut region_start = frame;
    let mut i = frame;

    while i + 3 <= seq.len() {
        let codon = &seq[i..i + 3];

        if stop_codons.contains(&codon) {
            let region_length = (i - region_start)/3;
            max_length = max_length.max(region_length);
            region_start = i + 3;
        }

        i += 3;
    }

    // Final region (after last stop, or entire frame if no stops)
    let final_region_length = (i - region_start) / 3;
    max_length.max(final_region_length)
}

/// Calculate the longest stop-free run in a given sequence, in codons
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

#[pyfunction]
fn calculate_stop_free_run_with_ratio(sequences: Vec<(String,String)>) -> PyResult<Vec<(String,usize, f64)>> {
    Ok(sequences
        .par_iter()
        .map(|(id, seq)| {
            let length = calculate_stop_free_run_single(seq);
            let ratio = if seq.is_empty() { 0.0 } else { length as f64 / seq.len() as f64 };
            (id.clone(), length, ratio)
        })
        .collect())
}

#[pyfunction]
fn calculate_gc_content(sequences: Vec<(String, String)>) -> PyResult<Vec<(String, f64)>> {
    Ok(sequences
        .par_iter()
        .map(|(id, seq)| {
            let gc_count = seq.chars().filter(|&c| c == 'G' || c == 'C' || c == 'g' || c == 'c').count();
            let gc_content = if seq.is_empty() { 0.0 } else { gc_count as f64 / seq.len() as f64 };
            (id.clone(), gc_content)
        })
        .collect())
}

#[pyfunction]
fn calculate_run_probability(run_lengths: Vec<(String, i32)>, gc_contents: Vec<(String, f64)>) -> PyResult<Vec<(String,  f64)>> {
    let gc_map: std::collections::HashMap<String, f64> = gc_contents.into_iter().collect();
    
    Ok(run_lengths
        .par_iter()
        .map(|(id, length)| {
            let gc_content = gc_map.get(id).cloned().unwrap_or(0.0);
            let probability = (1.0 - (((1.0 - gc_content).powi(2)*(1.0 + gc_content))/ 8.0) ).powi(*length); 
            (id.clone(), probability)
        })
        .collect())
}


/// A Rust implementation of stopFree
#[pymodule]
fn stopfree(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calculate_stop_free_run, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_stop_free_runs_with_ids, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_stop_free_run_with_ratio, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_run_probability, m)?)?;
    Ok(())
}



#[cfg(test)]
mod tests {
    use super::*;

    // ==================== Frame-level tests ====================

    #[test]
    fn test_frame0_no_stops() {
        let seq = "ATGATGATG"; // 9 nt = 3 codons
        assert_eq!(longest_stop_free_in_frame(seq, 0), 3);
    }

    #[test]
    fn test_frame0_stop_at_start() {
        let seq = "TAAATGATG"; // 9 nt
        // TAA|ATG|ATG → 0, then 6 = 2 codons
        assert_eq!(longest_stop_free_in_frame(seq, 0), 2);
    }

    #[test]
    fn test_frame0_stop_at_end() {
        let seq = "ATGATGTAA"; // 9 nt
        // ATG|ATG|TAA → 6, then 0 = 2 codons
        assert_eq!(longest_stop_free_in_frame(seq, 0), 2);
    }

    #[test]
    fn test_frame0_stop_in_middle() {
        let seq = "ATGTAAATG"; // 9 nt
        // ATG|TAA|ATG → 3, 3 = 1 codon
        assert_eq!(longest_stop_free_in_frame(seq, 0), 1);
    }

    #[test]
    fn test_frame0_multiple_stops() {
        let seq = "TAATAATAA"; // 9 nt
        // TAA|TAA|TAA → 0, 0, 0
        assert_eq!(longest_stop_free_in_frame(seq, 0), 0);
    }

    #[test]
    fn test_frame0_longest_region_at_end_exact_multiple() {
        // This is the Python bug case - at the frame level
        let seq = "TAGATGAAA"; // 9 nt
        // TAG|ATG|AAA → 0, then 6 = 2 codons
        assert_eq!(longest_stop_free_in_frame(seq, 0), 2);
    }

    #[test]
    fn test_frame0_longest_between_two_stops() {
        let seq = "TAAATGATGATGTAG"; // 15 nt
        // TAA|ATG|ATG|ATG|TAG → 0, 9, 0 = 3 codons
        assert_eq!(longest_stop_free_in_frame(seq, 0), 3);
    }

    #[test]
    fn test_frame1_finds_stop() {
        let seq = "ATAAGCC"; // 7 nt
        // Frame 1: TAA|GCC → 0, then 3 = 1 codon
        assert_eq!(longest_stop_free_in_frame(seq, 1), 1);
    }

    #[test]
    fn test_frame1_no_stop() {
        let seq = "AATGATG"; // 7 nt
        // Frame 1: ATG|ATG → 6 = 2 codons
        assert_eq!(longest_stop_free_in_frame(seq, 1), 2);
    }

    #[test]
    fn test_frame2_finds_stop() {
        let seq = "CCTAAGG"; // 7 nt
        // Frame 2: TAA|GG → 0, then 2
        assert_eq!(longest_stop_free_in_frame(seq, 2), 0);
    }

    #[test]
    fn test_frame_with_partial_codon_at_end() {
        let seq = "ATGATGA"; // 7 nt
        // Frame 0: ATG|ATG|A → 6 (partial A doesn't count) = 2 codons
        assert_eq!(longest_stop_free_in_frame(seq, 0), 2);
        // Frame 1: TGA|TGA → 0, 0
        assert_eq!(longest_stop_free_in_frame(seq, 1), 0);
        // Frame 2: GAT|GA → 3 = 1 codon
        assert_eq!(longest_stop_free_in_frame(seq, 2), 1);
    }

    #[test]
    fn test_frame_empty_sequence() {
        assert_eq!(longest_stop_free_in_frame("", 0), 0);
        assert_eq!(longest_stop_free_in_frame("", 1), 0);
        assert_eq!(longest_stop_free_in_frame("", 2), 0);
    }

    #[test]
    fn test_frame_too_short_for_codon() {
        assert_eq!(longest_stop_free_in_frame("AT", 0), 0);
        assert_eq!(longest_stop_free_in_frame("AT", 1), 0);
        assert_eq!(longest_stop_free_in_frame("AT", 2), 0);
    }

    #[test]
    fn test_all_stop_types_recognized_frame0() {
        assert_eq!(longest_stop_free_in_frame("TAA", 0), 0);
        assert_eq!(longest_stop_free_in_frame("TAG", 0), 0);
        assert_eq!(longest_stop_free_in_frame("TGA", 0), 0);
    }

    // ==================== Sequence-level tests (3 frames) ====================

    #[test]
    fn test_sequence_picks_best_frame() {
        let seq = "ATAAGCC"; // 7 nt
        // Frame 0: ATA|AGC|C → 6
        // Frame 1: TAA|GCC → 0, 3
        // Frame 2: AAG|CC → 3
        // Max = 6 (frame 0) = 2 codons
        assert_eq!(longest_stop_free_in_sequence(seq), 2);
    }

    #[test]
    fn test_sequence_all_frames_have_stops() {
        let seq = "TAATAGTGA"; // 9 nt
        // Frame 0: TAA|TAG|TGA → 0
        // Frame 1: AAT|AGT|GA → 6
        // Frame 2: ATA|GTG|A → 6
        // Max = 6 = 2 codons
        assert_eq!(longest_stop_free_in_sequence(seq), 2);
    }

    #[test]
    fn test_sequence_frame1_wins() {
        // Construct where frame 1 has longest run
        let seq = "TAAATGATGATG"; // 12 nt
        // Frame 0: TAA|ATG|ATG|ATG → 0, 9
        // Frame 1: AAA|TGA|TGA|TG → 3, 0, 0
        // Frame 2: AAT|GAT|GAT|G → 9
        // Max = 9 = 3 codons
        assert_eq!(longest_stop_free_in_sequence(seq), 3);
    }

    // ==================== Reverse complement tests ====================

    #[test]
    fn test_rc_basic() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
    }

    #[test]
    fn test_rc_stops_become_nonstops() {
        assert_eq!(reverse_complement("TAA"), "TTA");
        assert_eq!(reverse_complement("TAG"), "CTA");
        assert_eq!(reverse_complement("TGA"), "TCA");
    }

    #[test]
    fn test_rc_lowercase() {
        assert_eq!(reverse_complement("atgc"), "GCAT");
    }

    #[test]
    fn test_rc_preserves_n() {
        assert_eq!(reverse_complement("ANG"), "CNT");
    }

    // ==================== Full 6-frame tests ====================

    #[test]
    fn test_full_no_stops_anywhere() {
        let seq = "ATGATGATG";//= 3 codons
        // No stops in fwd or RC
        assert_eq!(calculate_stop_free_run_single(seq), 3);
    }

    #[test]
    fn test_full_stops_in_forward_rc_clean() {
        let seq = "TAATAGTGA"; // All stops in frame 0
        // Forward max = 6 (frames 1 or 2)
        // RC = TCACTATTA, no stops
        // RC max = 9 = 3 codons
        assert_eq!(calculate_stop_free_run_single(seq), 3);
    }

    #[test]
    fn test_full_forward_clean_stops_in_rc() {
        // TTA, CTA, TCA reverse complement TO stops
        let seq = "TTACTATTA"; // 9 nt
        // Forward: no stops → 9
        // RC = TAATAGTAA, has stops
        // Forward wins with 9 = 3 codons
        assert_eq!(calculate_stop_free_run_single(seq), 3);
    }

    #[test]
    fn test_full_palindrome_stops_both_directions() {
        let seq = "TAATTA"; // 6 nt, RC = TAATTA
        // Forward Frame 0: TAA|TTA → 0, 3
        // RC is identical
        // Max = 3 = 1 codon
        assert_eq!(calculate_stop_free_run_single(seq), 1);
    }

    #[test]
    fn test_full_case_insensitive() {
        assert_eq!(
            calculate_stop_free_run_single("taatagtga"),
            calculate_stop_free_run_single("TAATAGTGA")
        );
    }

    // ==================== Batch tests ====================

    #[test]
    fn test_batch() {
        let seqs = vec![
            "ATG".to_string(),
            "TAATAGTGA".to_string(),
        ];
        let results = calculate_stop_free_run(seqs).unwrap();
        assert_eq!(results, vec![1, 3]);
    }

    #[test]
    fn test_batch_with_ids() {
        let tuples = vec![
            ("s1".to_string(), "ATG".to_string()),
            ("s2".to_string(), "TAATAGTGA".to_string()),
        ];
        let results = calculate_stop_free_runs_with_ids(tuples).unwrap();
        assert_eq!(results[0], ("s1".to_string(), 1));
        assert_eq!(results[1], ("s2".to_string(), 3));
    }

    #[test]
    fn test_calculate_gc_content_50pc() {
        let tuples = vec![
            ("s1".to_string(), "ATGC".to_string()),
        ];
        let results = calculate_gc_content(tuples).unwrap();
        assert_eq!(results[0], ("s1".to_string(), 0.5));
    }

    #[test]
    fn test_calculate_gc_content_0pc() {
        let tuples = vec![
            ("s1".to_string(), "ATAT".to_string()),
        ];
        let results = calculate_gc_content(tuples).unwrap();
        assert_eq!(results[0], ("s1".to_string(), 0.0));
    }

    #[test]
    fn test_calculate_run_probability_k60_gc70pc(){
        // let tuple = vec![("s1".to_string(), "ATAGCGCGCG".to_string())];
        let run_lengths = vec![("s1".to_string(), 62)];
        let gc_contents = vec![("s1".to_string(), 0.50)];
        let results = calculate_run_probability(run_lengths, gc_contents).unwrap();
        println!("{:?}", results);
        println!("{}", (results[0].1 - 0.05096727139795341).abs());

        // Have to check for difference form machine epsilon for 64-bit floating point
        assert!((results[0].1 - 0.05096727139795341).abs() < 1.11e-16);
    }

}