# Pepars - Protein Engineering via PARallel Sequencing
Pepars is a Python package containing various utilities for dealing with parallel sequencing data (e.g. NGS FASTQ files) for protein engineering contexts.

## Installation
You can install Pepars via pip:
```
pip install git+https://github.com/GradinaruLab/pepars.git
```

## Functionality

Pepars is broken up into packages roughly based on functionality. The main
packages are:
- alignment
- analysis
- fileio
- plotting
- utils
- simulation

Details about the functionality of each package are below.

### Alignment

The alignment package contains functions for taking FASTQ files and extracting
variant regions. The main element of the alignment package is the Aligner class,
which is an abstract class that defines how alignment should function.

To perform alignment on a set of FASTQ files, instantiate one of the subclasses
of Aligner (e.g. Perfect_Match_Aligner or Bowtie_Aligner), and then call the
align function:

```
align(
    template,
    FASTQ_file_sets,
    alignment_parameters=None,
    progress_callback=None,
    reverse_complement_template=None):
```

See ```examples/alignment.ipynb``` for a typical use case.

The alignment parameters vary by aligner, as below:

#### Perfect_Match_Aligner Parameters
```
variant_sequence_quality_threshold=0,
mismatch_quality_threshold=0
```

#### Bowtie_Aligner Parameters
```
working_directory=os.getcwd(),
is_local=False,
output_frequency=1e5,
approach=None,
allow_insertions_deletions=False,
quality_threshold=0.0
```

#### Subclassing Aligner

To subclass the Aligner class, your class needs to implement the internal
_align method.

### Analysis

The analysis package has a variety of analyses relevant to parallel
sequencing-based protein engineering experiments. Some examples are:
- ```amino_acids.get_amino_acid_codon_biases```: Get the expected bias of amino
acids for a degenerate nucleotide sequence.
- ```sequencing_reads.get_nucleotide_distribution```: Get the distribution of
nucleotides in a FASTQ file
- ```confidence.get_sequence_confidences```: Get the normalized confidences of
a list of sequences and their counts, based on a few different confidence
metrics 

### Fileio

Some simple file wrappers, designed to make reading/writing CSV and sequence
count files easier.

### Plotting

A set of wrappers around Plotly, designed to make it take just one or a few
lines of code to get nice plots, either interactively or exported. Most of the
plotting functions in here internally call ```plotting.generate_plotly_plot```,
which is a catch-all Plotly wrapper that prints to screen, or writes to a file,
or both. To have plots generate interactively, make sure start your notebook
with ```plotting.init_notebook_mode()```

Some useful plots:

- ```plotting.plot_histogram```: A simple Plotly histogram wrapper
- ```plotting.plot_scatter```: A simple Plotly scatter plot wrapper
- ```plotting.plot_count_distribution```: Plot the count distribution of variant
sequences over multiple samples
- ```DNA.plot_amino_acid_bias```: Plot a heatmap of amino acid bias, given
sequence counts and a template

### Utils

A variety of potentially useful utility functions. Some highlights:

- ```DNA```: All the typical bioinformatics stuff: IUPAC grammar, the genetic
code, translation/complement functions
- ```FASTQ_File```: An easy-to-use FASTQ file iterator. Operates seemlessly on
both gzipped and raw FASTQ files, and lets you iterate by sequences, quality
scores, or both
- ```FASTQ_File_Set```: A convenient way to iterate line-by-line along multiple
FASTQ files in parallel - e.g. for paired end reads
- ```Sequence_Trie```: A fast storage and query data structure for sequence
information. Takes advantage of a trie structure to rapidly search for/iterate
over one-off or two-off mutants.
- ```AminoAcid```: A class for each Amino acid - useful for extracting physical
properties

### Simulation 

Currently, this is just some functions for generating random mutants and their
counts based on a template. Useful for testing.

## License
License information can be found in the LICENSE file