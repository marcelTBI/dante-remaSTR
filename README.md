# Dante (Remastered)

**Dante** ("Da Amazing NucleoTide Exposer") is an algorithm designed for genotyping STR alleles based on NGS reads originating from the STR locus of interest. This method takes into account natural deviations from the expected sequence, including variations in repeat count, sequencing errors, ambiguous bases, and complex loci containing various motifs.

The reported figures provide evidence for expanded alleles which are too long to be captured by a single NGS read, as well as allelic single point mutations, small insertions, and deletions that may be relevant for diagnostic evaluations.

## Project Structure

The entire algorithm is divided into two components:

1. **Annotation of Mapped Reads** - Using the `remastr` tool written in Rust for enhanced parallelism and speed.
   - Managed in a separate repository.
   - **Input:** BAM file with mapped reads, TXT file with HGVS nomenclature, FASTA file with reference.
   - **Output:** TSV file with annotated reads. Columns include:
     - **motif:** Motif nomenclature
     - **read_sn:** Number of read for each motif
     - **read_id:** Read identifier as found in BAM
     - **mate_order:** Information about read pairing: 1 - left read, 2 - right read, 0 - unpaired
     - **read:** Read sequence
     - **reference:** Aligned annotated motif in the read
     - **modules:** Aligned numbers of motif modules
     - **log_likelihood:** Log likelihood of the annotation

2. **Genotyping of Annotated Reads, Visualization, and Reporting** - Utilizes the `dante_remastr` tool in Python.
   - **Input:** TSV file with annotated reads (output from the `remastr` tool).
   - **Output:** TSV file with genotyped read motifs and confidences:
     - **motif_name:** Motif nomenclature
     - **motif_sequence:** Motif STR (part of the motif)
     - **chromosome:** Chromosome of the motif
     - **start:** Start location of the motif
     - **end:** End location of the motif
     - **allele1:** Number of repeats of the STR (or B - background, E - expanded) on the first allele
     - **allele2:** Number of repeats of the STR (or B - background, E - expanded) on the second allele
     - **confidence:** Prediction confidence in percent (calculated as a proportion of the likelihood of the predicted state versus all states)
     - **conf_allele1:** Confidence of the prediction for the first allele
     - **conf_allele2:** Confidence of the prediction for the second allele
     - **quality_reads:** Number of quality reads for this STR with both primers
     - **one_primer_reads:** Number of quality reads for this STR with exactly one primer
     - **filtered_reads:** Filtered reads mapped on this location (due to no primer, non-consecutive modules, too many errors, or insufficient repetitions or bases in modules)
     - **conf_background:** Confidence of the background prediction for both alleles
     - **conf_background_all:** Confidence of the background prediction for at least one allele
     - **conf_extended:** Confidence of the expanded allele prediction for both alleles
     - **conf_extended_all:** Confidence of the expanded allele prediction for at least one allele

## Getting Started

These instructions will help you set up the project on your local machine.

### Acquisition of Repository

Clone this repository using:

```bash
git clone https://github.com/marcelTBI/dante-remaSTR.git
```

### Conda Environment

Dante relies on libraries specified in the `conda_env.yaml`. Create the conda environment as follows (though `conda` works, `mamba` is recommended for speed):

```bash
mamba env create -n env_dante -f conda_env.yaml
```

### Download Remastr Subtool

```bash
git clone https://gitlab.com/andrejbalaz/remastr.git
```

#### Compilation and Running

The `remastr` Rust component needs to be compiled. Navigate to the submodule and run the compilation:

```bash
cd remastr/dante_cli
conda activate env_dante
cargo build --release  # Rust needs to be installed (https://www.rust-lang.org/tools/install)
```

Upon successful compilation, the release version of the `dante_cli` tool should appear in the `remastr/target/release` directory.

The `dante_remastr` Python part does not require compilation.

Run both parts as shown:

```bash
./remastr/target/release/dante_cli -b <reads_bam> -m <STRSet_file.tsv> -o remastr_output.tsv
python ./dante-remaSTR/dante_remastr_standalone.py -i remastr_output.tsv -o remastr_output -v
```

### Outputs

The tool can generate comprehensive visualizations and user reports in `.html` format using the `--verbose` or `-v` switch. However, this is feasible for up to about 100 motifs.

Notable result files include:

- `report.html` - Visual report in `.html` (generated only with the `-v` flag)
- `variants.tsv` - STRs in custom table format
- `variants.vcf` - STR variants in standard VCF format
- Folder `alignments` - Alignments of annotated reads per motif

### Example

An example is provided in the `example` directory, which can be run with:

```bash
cd example
conda activate env_dante
<path_to_remastr>/target/release/dante_cli -b father_sorted.bam -m STR_set.tsv -o father_remastr_output.tsv
python <path_to_danteSTR>/dante_remastr_standalone.py -i father_remastr_output.tsv -o father -v
```

Upon successful execution (ignore `Warning: ['23', 'B'] ['24', '11'] 22|13 24|11 is inconsistent.` warnings), you should be able to view the report in `father/report.html`.
