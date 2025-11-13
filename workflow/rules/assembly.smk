"""
Assembly Module
===============

Assembly rules for virome contigs.

Workflow:
1. Merge overlapping read pairs (bbmerge)
2. Concatenate all samples for coassembly OR assemble individually
3. Assemble using MEGAHIT
4. Filter and validate contigs

Input: Cleaned paired-end reads (from QC module or user-provided)
Output: Assembled contigs + assembly stats
"""

# ================================================================================
# Helper Functions for Entry Point Logic
# ================================================================================

def get_reads_for_merging(wildcards):
    """
    Get cleaned reads for bbmerge based on pipeline entry point.

    Returns paths to R1/R2 files depending on where pipeline starts.
    """
    start_from = config["pipeline"].get("start_from", "raw_reads")

    if start_from == "raw_reads":
        # Use QC module outputs (rrna_removed = final clean reads)
        return {
            "r1": f"{OUTDIR}/rrna_removed/{wildcards.sample}_R1.fastq.gz",
            "r2": f"{OUTDIR}/rrna_removed/{wildcards.sample}_R2.fastq.gz"
        }
    elif start_from == "cleaned_reads":
        # Use user-provided cleaned reads
        cleaned_dir = config["pipeline"].get("cleaned_reads_dir")
        if not cleaned_dir:
            raise ValueError(
                "Pipeline start_from='cleaned_reads' but 'cleaned_reads_dir' not specified in config!"
            )
        return {
            "r1": f"{cleaned_dir}/{wildcards.sample}_R1.fastq.gz",
            "r2": f"{cleaned_dir}/{wildcards.sample}_R2.fastq.gz"
        }
    else:
        raise ValueError(
            f"Invalid start_from value: '{start_from}'. "
            f"Must be 'raw_reads' or 'cleaned_reads'"
        )


# ================================================================================
# Read Merging with BBMerge
# ================================================================================

rule bbmerge:
    """
    Merge overlapping read pairs using BBMerge

    For short-insert libraries (common in viral metagenomics), many read pairs
    overlap and can be merged into single longer reads. This improves assembly
    quality and simplifies downstream analysis.

    Outputs:
    - merged: Successfully merged reads (single-end)
    - unmerged_r1/r2: Read pairs that couldn't be merged (paired-end)
    - hist: Insert size histogram for QC
    """
    input:
        unpack(get_reads_for_merging)
    output:
        merged = f"{OUTDIR}/bbmerge/{{sample}}_merged.fastq.gz",
        unmerged1 = f"{OUTDIR}/bbmerge/{{sample}}_R1_unmerged.fastq.gz",
        unmerged2 = f"{OUTDIR}/bbmerge/{{sample}}_R2_unmerged.fastq.gz",
        hist = f"{OUTDIR}/bbmerge/{{sample}}_hist.txt"
    log:
        f"{OUTDIR}/logs/bbmerge/{{sample}}.log"
    threads: 8
    resources:
        mem_mb = 16000
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        # Validate input files
        [ -s {input.r1} ] || {{ echo "Error: R1 file missing or empty" >> {log}; exit 1; }}
        [ -s {input.r2} ] || {{ echo "Error: R2 file missing or empty" >> {log}; exit 1; }}

        # Run bbmerge
        bbmerge.sh \
            in1={input.r1} \
            in2={input.r2} \
            out={output.merged} \
            outu1={output.unmerged1} \
            outu2={output.unmerged2} \
            ihist={output.hist} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2>&1 | tee {log}

        # Verify outputs
        gzip -t {output.merged} || {{ echo "Error: Merged output corrupted" >> {log}; exit 1; }}
        gzip -t {output.unmerged1} || {{ echo "Error: Unmerged R1 corrupted" >> {log}; exit 1; }}
        gzip -t {output.unmerged2} || {{ echo "Error: Unmerged R2 corrupted" >> {log}; exit 1; }}
        """


# ================================================================================
# Concatenation for Coassembly
# ================================================================================

rule concatenate_merged_reads:
    """Concatenate all merged reads for coassembly"""
    input:
        expand(f"{OUTDIR}/bbmerge/{{sample}}_merged.fastq.gz", sample=SAMPLES)
    output:
        f"{OUTDIR}/assembly/all_merged.fastq.gz"
    log:
        f"{OUTDIR}/logs/assembly/concat_merged.log"
    threads: 4
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        echo "Concatenating {input}" > {log}
        zcat {input} | gzip -c > {output} 2>> {log}
        gzip -t {output} && echo "✅ Output verified" >> {log}
        """


rule concatenate_unmerged_r1:
    """Concatenate all unmerged R1 reads for coassembly"""
    input:
        expand(f"{OUTDIR}/bbmerge/{{sample}}_R1_unmerged.fastq.gz", sample=SAMPLES)
    output:
        f"{OUTDIR}/assembly/all_unmerged_R1.fastq.gz"
    log:
        f"{OUTDIR}/logs/assembly/concat_r1.log"
    threads: 4
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        echo "Concatenating {input}" > {log}
        zcat {input} | gzip -c > {output} 2>> {log}
        gzip -t {output} && echo "✅ Output verified" >> {log}
        """


rule concatenate_unmerged_r2:
    """Concatenate all unmerged R2 reads for coassembly"""
    input:
        expand(f"{OUTDIR}/bbmerge/{{sample}}_R2_unmerged.fastq.gz", sample=SAMPLES)
    output:
        f"{OUTDIR}/assembly/all_unmerged_R2.fastq.gz"
    log:
        f"{OUTDIR}/logs/assembly/concat_r2.log"
    threads: 4
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        echo "Concatenating {input}" > {log}
        zcat {input} | gzip -c > {output} 2>> {log}
        gzip -t {output} && echo "✅ Output verified" >> {log}
        """


# ================================================================================
# Assembly with MEGAHIT
# ================================================================================

rule megahit_coassembly:
    """
    Coassembly of all samples using MEGAHIT

    Combines all merged (single-end) and unmerged (paired-end) reads into a
    single assembly. This is recommended for virome studies to:
    - Increase coverage of low-abundance viruses
    - Enable strain-level resolution across samples
    - Reduce computational cost vs per-sample assembly

    Uses MEGAHIT's default k-mer list optimized for metagenomes.
    """
    input:
        merged = f"{OUTDIR}/assembly/all_merged.fastq.gz",
        r1 = f"{OUTDIR}/assembly/all_unmerged_R1.fastq.gz",
        r2 = f"{OUTDIR}/assembly/all_unmerged_R2.fastq.gz"
    output:
        contigs = f"{OUTDIR}/assembly/megahit/final.contigs.fa",
        done = touch(f"{OUTDIR}/assembly/megahit/.done")
    params:
        out_dir = f"{OUTDIR}/assembly/megahit",
        min_contig = config.get("assembly", {}).get("min_contig_length", 1000)
    log:
        f"{OUTDIR}/logs/assembly/megahit.log"
    threads: 24
    resources:
        mem_mb = 100000  # 100GB for large coassemblies
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        # Remove output directory if it exists (megahit won't overwrite)
        rm -rf {params.out_dir}

        # Run megahit
        megahit \
            -r {input.merged} \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.out_dir} \
            --min-contig-len {params.min_contig} \
            --k-list 21,29,39,59,79,99,119,141 \
            -t {threads} \
            &> {log}

        # Verify output
        [ -f {output.contigs} ] || {{ echo "Error: Assembly failed" >> {log}; exit 1; }}
        """


rule megahit_individual:
    """
    Individual sample assembly using MEGAHIT

    Assembles each sample separately. Use this when:
    - Samples are from different conditions/timepoints
    - You need sample-specific strain resolution
    - Coassembly is too memory-intensive
    """
    input:
        merged = f"{OUTDIR}/bbmerge/{{sample}}_merged.fastq.gz",
        r1 = f"{OUTDIR}/bbmerge/{{sample}}_R1_unmerged.fastq.gz",
        r2 = f"{OUTDIR}/bbmerge/{{sample}}_R2_unmerged.fastq.gz"
    output:
        contigs = f"{OUTDIR}/assembly/per_sample/{{sample}}/final.contigs.fa",
        done = touch(f"{OUTDIR}/assembly/per_sample/{{sample}}/.done")
    params:
        out_dir = f"{OUTDIR}/assembly/per_sample/{{sample}}",
        min_contig = config.get("assembly", {}).get("min_contig_length", 1000)
    log:
        f"{OUTDIR}/logs/assembly/megahit_{{sample}}.log"
    threads: 12
    resources:
        mem_mb = 50000  # 50GB per sample
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        # Remove output directory if it exists
        rm -rf {params.out_dir}

        # Run megahit
        megahit \
            -r {input.merged} \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.out_dir} \
            --min-contig-len {params.min_contig} \
            --k-list 21,29,39,59,79,99,119,141 \
            -t {threads} \
            &> {log}

        # Verify output
        [ -f {output.contigs} ] || {{ echo "Error: Assembly failed" >> {log}; exit 1; }}
        """


# ================================================================================
# Assembly Statistics and QC
# ================================================================================

rule assembly_stats:
    """
    Calculate assembly statistics

    Metrics:
    - Number of contigs
    - Total assembly size
    - N50, L50
    - Longest contig
    - GC content
    """
    input:
        contigs = f"{OUTDIR}/assembly/megahit/final.contigs.fa" if config["pipeline"].get("assembly_strategy", "coassembly") == "coassembly"
                  else expand(f"{OUTDIR}/assembly/per_sample/{{sample}}/final.contigs.fa", sample=SAMPLES)
    output:
        stats = f"{OUTDIR}/reports/assembly_stats.tsv"
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/assembly_stats.py"
