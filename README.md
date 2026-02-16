# stopfree

[![Rust Tests](https://github.com/afg1/stopfree/actions/workflows/cargo-test.yml/badge.svg)](https://github.com/afg1/stopfree/actions/workflows/cargo-test.yml)

This is a rust implementation of the stopFree coding potential tool defined in [Evaluating computational tools for protein-coding sequence detection: Are they up to the task?](https://www.biorxiv.org/content/10.1101/2024.05.16.594598v3.full.pdf)

It measures the longest stop codon free length of sequence in each of the 6 possible reading frames.

## Installation

You should be able to just 

```
pip install stopfree
```

because there's a bunch of pre-built wheels available on pypi.

## Usage

To calculate stop codon free lengths, you need to build a list of tuples of ID and sequence, like this:

```
all_fasta = []
for idx, item in enumerate(read_fasta("mgi.fasta")):
    seqid = item.defline
    seq = item.sequence
    all_fasta.append((seqid, seq))
```
(N.B. Here I used the `fasta_reader` package to read a FASTA file into an iterable)

Then call the stopfree functions like this:

```
stop_free_run_lengths = stopfree.calculate_stop_free_runs_with_ids(all_fasta) ## Gives the number of codons with no stop codon in
gc_contents = stopfree.calculate_gc_content(all_fasta) ## Calculates the per-sequence GC content for each sequence (not across all!)
run_probabilities = stopfree.calculate_run_probability(stop_free_run_lengths, gc_contents) ## Calculated the probability of a run of the observed length
```

### Probability calculation

We use an iid assumption, with correction for GC content of the sequence. With this, the probability of a stop codon is given by

```
p_stop = ((1 - p_gc)**2 * (1 + p_gc))/8
```
where `p_gc` is the gc content fraction.

So the probability of observing no stop codon in a run of k codons is therefore:

```
p_nostop = (1 - p_stop)**k
```

You can set whatever threshold for 'protein coding' you like then, based on this probability.



### Better probability

The iid assumption doesn't quite account for everything. To do a really stellar job, we would want to include dinucleotide biases, and the fact that the 6 reading frames are correlated. This is some future work though for now.

