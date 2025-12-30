# Lab Virome QC Pipeline

A comprehensive quality control and assembly pipeline for VLP-enriched virome sequencing data generated from RdAB (Random displacement Amplification) protocol and Illumina NovaSeq sequencing.

---

## Overview

This Snakemake pipeline provides robust QC specifically designed for:

- **VLP (Virus-Like Particle) enriched samples**
- **RdAB amplification** (RT + Random priming + PCR)
- **Illumina NovaSeq sequencing** (2-channel chemistry)
- **Mechanical shearing or tagmentation-based** library prep

### Key Features

✅ **NovaSeq-specific QC** - PolyG tail removal (critical for 2-channel chemistry)
✅ **VLP enrichment assessment** - ViromeQC enrichment scoring
✅ **Contamination flagging & removal** - PhiX/vector flagging (non-destructive), host & rRNA removal
✅ **Statistical outlier detection** - IQR-based contamination QC with publication-quality visualizations
✅ **Optical duplicate removal** - Illumina patterned flow cell artifacts
✅ **Automated QC flagging** - Pass/fail criteria for each sample
✅ **Rich reporting** - MultiQC dashboard with all metrics
✅ **Modular assembly** - Optional viral metagenome assembly (individual or coassembly strategies)
✅ **Flexible entry points** - Start from raw reads or previously cleaned reads

---

## Pipeline Workflow

### Full Pipeline (QC + Assembly)

```
Raw Reads (NovaSeq FASTQ)
    ↓
[1] FastQC (raw reads)
    ↓
[2] Clumpify (remove optical duplicates)
    ↓
[3] fastp (CRITICAL: polyG removal + adapter trimming + QC)
    ↓
[4] FastQC (trimmed reads)
    ↓
[5] BBDuk (PhiX/vector contamination flagging) ← NEW: Non-destructive detection
    ↓
[6] minimap2 (host depletion) ← QC metric for VLP success
    ↓
[7] ViromeQC (enrichment assessment) ← PRIMARY QC metric
    ↓
[8] BBDuk (rRNA removal)
    ↓
[9] FastQC (final clean reads)
    ↓
[10] MultiQC (aggregate all reports)
    ↓
[11] Contamination analysis (statistical outlier detection + plots)
    ↓
Clean reads + QC reports + Sample flags + Contamination plots
    ↓
    ↓ [OPTIONAL: Assembly Module]
    ↓
[12] BBMerge (merge overlapping read pairs)
    ↓
[13] MEGAHIT assembly
    │   ├─→ Individual assembly (per-sample assemblies)
    │   └─→ Coassembly (all samples pooled)
    ↓
Assembled contigs + Assembly statistics
```

**Note:** The pipeline is modular - you can run QC-only or QC + Assembly depending on your needs.

---

## Installation

### Requirements

- [Conda/Mamba](https://github.com/conda-forge/miniforge) (for environment management)
- [Snakemake](https://snakemake.readthedocs.io/) ≥7.0

### Quick Start

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/lab-virome-QC.git
cd lab-virome-QC

# Install Snakemake (if not already installed)
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake

# Test installation
snakemake --version
```

---

## Setup

### 1. Prepare Reference Databases

**rRNA Database Requirements**

The rRNA removal step requires a SILVA database that includes both:
- **SSU rRNA** (16S/18S - small subunit ribosomal RNA)
- **LSU rRNA** (23S/28S - large subunit ribosomal RNA)

Note: SSU-only databases (commonly used for 16S amplicon sequencing) will not remove LSU rRNA contamination in virome samples.

**Manual download and combine SILVA databases:**

```bash
# Download SILVA 138.1 SSU (16S/18S) - NR99 version
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

# Download SILVA 138.1 LSU (23S/28S) - NR99 version
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz

# Combine into one comprehensive database
zcat SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz > resources/silva_rrna.fasta
```

---

**Quick Setup (Recommended):**

Download all required reference databases with one command:

```bash
bash scripts/setup_references.sh human
```

This will download:
- PhiX174 reference (~5 KB)
- Human genome GRCh38 (~900 MB compressed, ~3 GB uncompressed)
- SILVA rRNA database (~150 MB compressed, ~350 MB uncompressed)

**Time:** 15-30 minutes depending on internet speed
**Disk space:** ~4-5 GB total

**For other organisms:**
```bash
bash scripts/setup_references.sh mouse    # For mouse samples
bash scripts/setup_references.sh custom   # For custom genome (will prompt for URL)
```

**Manual/Individual Downloads:**

If you prefer to download databases individually:

```bash
bash resources/download_phix.sh           # PhiX174 reference
bash resources/download_host.sh human     # Host genome
bash resources/download_silva.sh          # SILVA rRNA database
```

**For detailed documentation, troubleshooting, and advanced options:**
See `resources/README.md`

### 2. Configure Sample Information

You have three options for specifying samples:

**Option A: Auto-detect samples from directories (recommended for many samples)**

The pipeline supports flexible auto-detection with three modes:

**A1. Single directory (simple case)**
```yaml
sample_auto_detection:
  enabled: true
  input_dir: "data/raw"              # Directory containing your FASTQ files
  r1_pattern: "*_R1.fastq.gz"        # Pattern for R1 files
  r2_pattern: "*_R2.fastq.gz"        # Pattern for R2 files
```

**A2. Multiple directories (files spread across locations) - NEW!**
```yaml
sample_auto_detection:
  enabled: true
  input_dirs:                        # Multiple directories
    - "data/raw"
    - "data/additional_samples"
    - "/path/to/external/data"
  r1_pattern: "*_R1.fastq.gz"
  r2_pattern: "*_R2.fastq.gz"
  conflict_resolution: "prefix_dir"  # Handle duplicate sample names
```

**A3. Recursive scanning (files in subdirectories) - NEW!**
```yaml
sample_auto_detection:
  enabled: true
  input_dir: "data"                  # Root directory to scan
  recursive: true                    # Scan all subdirectories
  max_depth: 3                       # Optional: limit recursion depth
  r1_pattern: "*_R1.fastq.gz"
  r2_pattern: "*_R2.fastq.gz"
  conflict_resolution: "prefix_dir"  # Recommended for recursive mode
```

**Conflict resolution strategies** (for duplicate sample names):
- `"error"` - Stop with error if duplicates found (default)
- `"prefix_dir"` - Add directory name prefix (e.g., `dir1_sample1`, `dir2_sample1`)
- `"newest"` - Use sample from most recently modified directory

Common file naming patterns:
- `*_R1.fastq.gz` / `*_R2.fastq.gz` (default)
- `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz` (Illumina default naming)
- `*_1.fq.gz` / `*_2.fq.gz` (short form)
- `*.R1.fastq.gz` / `*.R2.fastq.gz` (dot separator)

**Example configurations:**
- Single directory: `config/config_auto_detect_example.yaml`
- Multiple directories: `config/config_multi_directory_example.yaml`
- Recursive scanning: `config/config_recursive_example.yaml`

**Option B: Edit config.yaml directly**
```yaml
samples:
  sample1:
    r1: "data/raw/sample1_R1.fastq.gz"
    r2: "data/raw/sample1_R2.fastq.gz"
  sample2:
    r1: "data/raw/sample2_R1.fastq.gz"
    r2: "data/raw/sample2_R2.fastq.gz"
```

**Option C: Use sample sheet**
```bash
# Edit config/samples.tsv
sample	r1	r2
sample1	data/raw/sample1_R1.fastq.gz	data/raw/sample1_R2.fastq.gz
sample2	data/raw/sample2_R1.fastq.gz	data/raw/sample2_R2.fastq.gz
```

### 3. Adjust QC Thresholds (Optional)

Edit `config/config.yaml` to set custom QC thresholds based on your lab's historical data:

```yaml
qc_thresholds:
  min_enrichment_score: 10      # ViromeQC enrichment score (adjust for your preps)
  max_host_percent: 10          # Maximum % host reads (adjust for your preps)
  max_rrna_percent: 20          # Maximum % rRNA after removal
  min_final_reads: 100000       # Minimum reads after QC (adjust for sequencing depth)
```

**Note:** These are example starting points. Adjust based on your lab's VLP preparation performance and sequencing platform.

### 4. Resource Requirements

The pipeline is designed for NovaSeq data and requires adequate computational resources:

**Memory Requirements:**

| Step | Memory | Notes |
|------|--------|-------|
| Optical duplicate removal | 48 GB | Scales with read count (~0.5GB per million reads) |
| Host depletion | 24 GB | Minimap2 index + read processing |
| rRNA removal | 32 GB | 5GB for SILVA SSU+LSU database + read processing |
| ViromeQC | 16 GB | Bowtie2/Diamond alignments |

These allocations handle NovaSeq samples up to 100 million reads. For smaller datasets or resource-limited systems, you can override these values using cluster profiles (see `profile/slurm/` for examples).

**Computational Resources:**
- **Threads:** Most rules use 4-8 threads
- **Runtime:** ~2-4 hours per sample (depends on sample size and cluster load)
- **Storage:** ~10-20 GB per sample for intermediate files

**Assembly Resource Requirements (if enabled):**

| Step | Memory | Threads | Notes |
|------|--------|---------|-------|
| BBMerge | 16 GB | 8 | Merge overlapping read pairs |
| MEGAHIT (individual) | 16-32 GB | 12-16 | Per-sample assembly |
| MEGAHIT (coassembly) | 64-128 GB | 24 | Memory scales with total data |

---

### 5. Pipeline Modularity and Assembly Configuration

The pipeline supports flexible entry points and optional assembly:

#### Pipeline Modes

**Mode 1: QC Only (default)**
```yaml
pipeline:
  run_assembly: false  # or omit this line
```
Runs quality control only, outputs clean reads.

**Mode 2: QC + Assembly**
```yaml
pipeline:
  run_assembly: true
  assembly_strategy: "individual"  # or "coassembly"
```
Runs complete QC followed by viral metagenome assembly.

**Mode 3: Assembly from Existing Clean Reads**
```yaml
pipeline:
  start_from: "cleaned_reads"
  cleaned_reads_dir: "/path/to/clean_reads"
  run_assembly: true
  assembly_strategy: "individual"
```
Skips QC, uses existing clean reads for assembly. Useful for re-assembly with different parameters or trying different assembly strategies.

#### Assembly Strategy Options

**Individual Assembly** (recommended for most virome studies)
```yaml
assembly_strategy: "individual"
```

**Advantages:**
- Preserves sample-specific viral variants
- Better resolution of dominant strains
- Captures sample-specific diversity
- Faster per-sample processing (can parallelize)

**When to use:**
- Comparing viral populations across conditions/timepoints
- Low sample counts (<10 samples)
- Samples with very different viral communities
- Need sample-specific variant information

**Coassembly** (for shared viral populations)
```yaml
assembly_strategy: "coassembly"
```

**Advantages:**
- Better assembly of shared/common viruses
- Higher effective coverage for low-abundance viruses
- Merges redundant sequences across samples
- Better for low-depth samples

**When to use:**
- Technical replicates or similar samples
- Longitudinal samples from same individual
- Batch assembly of related samples
- Reference catalog creation

#### Example Configurations

**Example 1: Full QC + Individual Assembly**
```yaml
# config/config.yaml
pipeline:
  run_assembly: true
  assembly_strategy: "individual"

samples:
  sample1:
    r1: "data/raw/sample1_R1.fastq.gz"
    r2: "data/raw/sample1_R2.fastq.gz"
```

**Example 2: Coassembly from Pre-QC'd Reads**
```yaml
# config/config.yaml
pipeline:
  start_from: "cleaned_reads"
  cleaned_reads_dir: "previous_run/clean_reads"
  run_assembly: true
  assembly_strategy: "coassembly"

sample_auto_detection:
  enabled: true
  input_dir: "previous_run/clean_reads"
  r1_pattern: "*_R1.fastq.gz"
  r2_pattern: "*_R2.fastq.gz"
```

---

## Usage

### Run Complete Pipeline

```bash
# Dry run (check what will be executed)
snakemake --use-conda -n

# Run locally with 8 cores
snakemake --use-conda --cores 8

# Run on HPC with SLURM (example: HTCF at WashU)
snakemake --profile profiles/htcf --jobs 50
```

**For HPC/cluster users:**
- WashU HTCF users: See [docs/HTCF_USAGE.md](docs/HTCF_USAGE.md) for complete setup
- Other clusters: The HTCF profile serves as a working example - see the [Adapting for Your Own Cluster](docs/HTCF_USAGE.md#adapting-for-your-own-cluster) section

**Quick HTCF Test:**
```bash
# On HTCF, activate environment and run test
source /ref/sahlab/software/miniforge3/bin/activate
conda activate snakemake_tutorial
bash scripts/test_htcf_run.sh
```

### Run Specific Steps

```bash
# Just run FastQC on raw reads
snakemake --use-conda --cores 4 results/fastqc/raw/

# Run up to adapter trimming
snakemake --use-conda --cores 8 results/fastp/

# Generate MultiQC report only
snakemake --use-conda --cores 2 results/multiqc/multiqc_report.html
```

### Generate DAG Visualization

```bash
snakemake --dag | dot -Tpng > dag.png
```

---

## Output Structure

### QC-Only Mode

```
results/
├── fastqc/                    # FastQC reports (raw, trimmed, final)
├── clumpify/                  # Optical duplicate removal
├── fastp/                     # Adapter trimming + QC
├── contamination_flagging/    # PhiX and vector contamination detection stats
│   ├── phix/                  # Per-sample PhiX contamination stats
│   └── univec/                # Per-sample vector/plasmid contamination stats
├── host_depleted/             # Host-depleted reads
├── rrna_removed/              # rRNA-depleted reads (clean)
├── viromeqc/                  # ViromeQC enrichment scores
├── clean_reads/               # Symlinks to final clean reads
├── reports/
│   ├── read_counts.tsv        # Read counts at each step
│   ├── sample_qc_flags.tsv    # Pass/fail flags per sample
│   ├── contamination_summary.tsv        # Contamination levels per sample
│   ├── contamination_bars.png           # Bar plot with outliers highlighted
│   ├── contamination_boxes.png          # Distribution box plots
│   ├── contamination_scatter.png        # PhiX vs vector correlation
│   └── contamination_heatmap.png        # Heatmap overview
├── multiqc/
│   └── multiqc_report.html    # Comprehensive QC dashboard
└── logs/                      # All log files
```

### With Assembly Enabled

```
results/
├── [All QC outputs above]
├── bbmerge/                   # Merged and unmerged read pairs
│   ├── {sample}_merged.fastq.gz        # Successfully merged reads
│   ├── {sample}_R1_unmerged.fastq.gz   # R1 reads that couldn't merge
│   ├── {sample}_R2_unmerged.fastq.gz   # R2 reads that couldn't merge
│   └── {sample}_hist.txt               # Insert size histogram
├── assembly/
│   ├── final.contigs.fa       # Final assembled contigs (coassembly)
│   │   OR
│   ├── per_sample/            # Individual assemblies (if assembly_strategy: "individual")
│   │   ├── {sample1}/
│   │   │   ├── final.contigs.fa         # Sample-specific contigs
│   │   │   └── intermediate_contigs/    # MEGAHIT k-mer intermediates
│   │   ├── {sample2}/
│   │   │   └── ...
│   └── megahit/               # MEGAHIT working directory (coassembly only)
└── reports/
    ├── assembly_stats.tsv     # Assembly statistics (N50, total size, etc.)
    └── [other QC reports]
```

---

## Interpreting Results

### 1. MultiQC Report

Open `results/multiqc/multiqc_report.html` in a web browser.

**Key sections to check:**
- **FastQC**: Look for polyG in overrepresented sequences (should be gone after fastp)
- **fastp**: Check adapter content, insert size distribution, duplication rates
- **Read counts**: Track reads retained at each step

### 2. Sample QC Flags

Check `results/reports/sample_qc_flags.tsv`:

```
sample    enrichment_score  pass_enrichment  pass_host  pass_rrna  overall_pass  notes
sample1   25.3              PASS             PASS       PASS       PASS          None
sample2   3.2               FAIL             FAIL       PASS       FAIL          Low_enrichment(3.2);High_host(15.3%)
```

**Interpretation:**
- **PASS**: Sample meets all QC criteria
- **FAIL**: Sample fails one or more criteria (check notes)
- **WARNING**: Insufficient data to assess

### 3. ViromeQC Enrichment Score

**Most important QC metric for VLP samples!**

- **Score >>1**: Good VLP enrichment (viral reads enriched vs bacterial)
- **Score ~1**: Poor enrichment (essentially a metagenome)
- **Score <1**: VLP prep failed (less viral content than typical metagenome)

Low scores indicate:
- Failed VLP preparation
- Excessive bacterial contamination
- Sample may need to be excluded or re-processed

Compare your samples to lab standards and previous successful runs to determine acceptable thresholds.

### 4. Contamination Flagging (PhiX & Vector)

**NEW in this version!** See detailed guide: [CONTAMINATION_FLAGGING.md](CONTAMINATION_FLAGGING.md)

Check contamination plots in `results/reports/`:
- **contamination_bars.png** - Sample-by-sample contamination levels (outliers in red)
- **contamination_boxes.png** - Distribution across your batch
- **contamination_scatter.png** - PhiX vs vector correlation
- **contamination_heatmap.png** - Overview with outliers boxed

**Key Points:**
- Uses **IQR-based outlier detection** (not fixed thresholds)
- Identifies samples that deviate from your batch median
- **Non-destructive**: Flags contamination but doesn't remove reads
- VLP samples typically have very low contamination (<0.01%)
- Outliers may indicate library prep issues (not VLP failure)

Check `results/reports/contamination_summary.tsv` for exact percentages.

### 5. Host Contamination

Check `results/host_depleted/*_host_stats.txt`

**Key concept:** High host reads indicate VLP prep failure (unlike metagenomes where host removal is routine cleanup).

Compare your samples to:
- Lab historical data for typical VLP prep performance
- Within-run median to identify outliers
- Manufacturer specifications if using commercial VLP enrichment kits

### 6. Read Retention

Check `results/reports/read_counts.tsv`

**Read retention varies widely** depending on:
- Library quality (adapter content, quality scores)
- NovaSeq polyG artifact prevalence (2-channel chemistry)
- rRNA contamination levels (especially for RT-based protocols)
- Sample type and VLP enrichment efficiency

Compare within-run samples to identify outliers with unusually high losses.

### 7. Assembly Results (if enabled)

#### BBMerge Statistics

Check `results/bbmerge/{sample}_hist.txt` for insert size distributions:

```
#Mean    206.9
#Median  213
#Mode    261
#STDev   52.5
#PercentOfPairs  55.241
```

**Key metrics:**
- **PercentOfPairs**: Merge rate (typical range: 45-75% for virome data)
  - Higher merge rates = shorter inserts = more overlapping reads
  - Lower merge rates = longer inserts = less overlap (both contribute to assembly)
- **Mean/Median insert size**: Should be 150-250 bp for 2x151bp sequencing
- **Mode**: Most common insert size in your library

**Interpretation:**
- Merge rates 50-70% are excellent for virome assembly
- Consistent insert sizes across samples indicate good library prep
- Wide standard deviation (>60 bp) may indicate heterogeneous fragmentation

#### Assembly Statistics

Check `results/reports/assembly_stats.tsv`:

```
sample      num_contigs  total_size_bp  mean_length_bp  longest_contig_bp  n50    l50  gc_percent
sample1     214          771425         3604            58285              6133   28   41.43
coassembly  1046         3132678        2994            93393              4598   121  43.50
```

**Key metrics:**

**N50** (most important assembly quality metric)
- Length-weighted median contig size
- Higher N50 = better assembly contiguity
- **Typical ranges for virome data:**
  - Individual assemblies: 2-85 kb (varies by sample complexity)
  - Coassembly: 3-10 kb (depends on shared viral content)

**Number of contigs**
- Individual: 30-500+ contigs (depends on viral diversity)
- Coassembly: 500-2000+ contigs (captures pan-virome)
- More contigs doesn't mean better/worse - reflects biological diversity

**Longest contig**
- Individual: 50-200 kb (often near-complete viral genomes)
- Coassembly: 50-150 kb
- Very long contigs (>100 kb) may represent:
  - Complete phage genomes (typical: 30-200 kb)
  - Large DNA viruses
  - Rare bacterial contamination (if GC% is ~50-65%)

**Total assembly size**
- Individual: 300 kb - 2 Mb per sample
- Coassembly: 2-10 Mb (depends on number of samples and diversity)
- Size reflects captured viral diversity, not necessarily sample quality

**GC content**
- Viromes typically: 35-50% GC
- High GC (>55%) may indicate bacterial contamination
- Very low GC (<30%) may indicate AT-rich viruses (e.g., some ssDNA viruses)

#### Comparing Assembly Strategies

Use your test data to decide between strategies:

**Individual assembly produced longer contigs (higher N50)?**
- Samples have dominant strain variants
- Use individual assembly for this dataset

**Coassembly produced longer contigs?**
- Samples share common viral populations
- Coassembly provides better resolution
- Consider coassembly for this dataset

**Similar results?**
- Either strategy is appropriate
- Choose based on downstream analysis needs

---

## Critical NovaSeq + VLP Considerations

### PolyG Tails (NovaSeq Artifact)

**Problem:** NovaSeq 2-channel chemistry calls dark cycles as high-quality G bases
**Solution:** fastp with `--trim_poly_g` (enabled by default in this pipeline)
**Check:** FastQC on raw reads should show GGGG in overrepresented sequences; should be gone after fastp

### VLP Enrichment Assessment

**Problem:** VLP prep can fail, resulting in metagenome-like contamination
**Solution:** ViromeQC quantifies enrichment score
**Action:** Compare scores across your batch to identify failed preps (significantly lower than batch median)

### Host Contamination as QC Metric

**Problem:** High host reads indicate VLP prep failure
**Solution:** Track % host reads; flag >10%
**Note:** Unlike metagenome QC, host contamination is QC failure, not just data cleanup

---

## Troubleshooting

### Pipeline fails at ViromeQC

**Error:** `viromeQC.py: command not found`

**Solution:**
```bash
# Install ViromeQC manually
conda activate viromeqc
pip install viromeqc
```

### High memory usage

**Problem:** BBTools rules consuming too much RAM

**Solution:** Reduce memory allocation in Snakefile:
```python
resources:
    mem_mb = 8000  # Reduce from 16000
```

### Pipeline is slow

**Problem:** Single-threaded execution

**Solution:**
```bash
# Use more cores
snakemake --use-conda --cores 16

# Or run on cluster
snakemake --use-conda --cluster "sbatch" --jobs 50
```

### Adapter sequences not removed

**Problem:** Custom primers not detected by fastp

**Solution:** Add custom adapter sequences to fastp rule:
```bash
--adapter_sequence CUSTOM_PRIMER_B_SEQ \
--adapter_sequence_r2 CUSTOM_PRIMER_B_R2_SEQ
```

---

## Citation

If you use this pipeline in your research, please cite:

- **Snakemake:** Mölder et al. (2021) F1000Research
- **fastp:** Chen et al. (2018) Bioinformatics
- **ViromeQC:** Zolfo et al. (2019) Nature Biotechnology
- **BBTools:** Bushnell (2014) BBMap
- **MultiQC:** Ewels et al. (2016) Bioinformatics

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contributing

We welcome contributions from lab members and the community!

### For Lab Members

- **New to the project?** Start with [LEARNING_RESOURCES.md](docs/LEARNING_RESOURCES.md)
- **Ready to contribute?** Read [CONTRIBUTING.md](CONTRIBUTING.md)
- **Need help with GitHub?** See [GITHUB_SETUP.md](docs/GITHUB_SETUP.md)
- **Questions?** Check our [Code of Conduct](CODE_OF_CONDUCT.md)

### Quick Start for Contributors

```bash
# Fork and clone
git clone https://github.com/YOUR_USERNAME/lab-virome-QC.git
cd lab-virome-QC

# Create a branch
git checkout -b feature/my-feature

# Make changes, commit, push
git add .
git commit -m "feat: description of changes"
git push origin feature/my-feature

# Open a Pull Request on GitHub
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

---

## Contact

For questions or issues:
- Open an [Issue](https://github.com/shandley/lab-virome-QC/issues)
- Email: scott.handley@wustl.edu
- Lab Slack: #virome-qc

---

## Acknowledgments

Developed for VLP-enriched virome analysis workflows.

Special considerations for:
- NovaSeq 2-channel chemistry artifacts
- RdAB amplification biases
- VLP enrichment quality assessment
- Mechanical shearing library preparation

---

## References

1. Chen et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884-i890.
2. Zolfo et al. (2019). Detecting contamination in viromes using ViromeQC. *Nature Biotechnology* 37:1408-1412.
3. Bushnell (2014). BBMap: A Fast, Accurate, Splice-Aware Aligner. *Lawrence Berkeley National Lab*
4. Ewels et al. (2016). MultiQC: summarize analysis results for multiple tools and samples. *Bioinformatics* 32(19):3047-3048.
5. 2024 benchmarking studies on virome QC best practices
