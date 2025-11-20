#!/usr/bin/env python3
"""
Calculate assembly statistics from FASTA contigs

Computes metrics including:
- Number of contigs
- Total assembly size (bp)
- N50, L50 (standard assembly quality metrics)
- Longest contig
- Mean contig length
- GC content (%)
"""

import pandas as pd
from pathlib import Path
from Bio import SeqIO
import numpy as np


def parse_fasta(fasta_file):
    """
    Parse FASTA file and return list of sequence records

    Args:
        fasta_file: Path to FASTA file

    Returns:
        list: List of SeqRecord objects
    """
    return list(SeqIO.parse(fasta_file, "fasta"))


def calculate_n50_l50(lengths):
    """
    Calculate N50 and L50 metrics

    N50: Length of the shortest contig for which longer and equal contigs
         cover at least 50% of the assembly
    L50: Number of contigs whose summed length is N50

    Args:
        lengths: List of contig lengths

    Returns:
        tuple: (N50, L50)
    """
    # Sort lengths in descending order
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    half_length = total_length / 2

    cumsum = 0
    n50 = 0
    l50 = 0

    for length in sorted_lengths:
        cumsum += length
        l50 += 1
        if cumsum >= half_length:
            n50 = length
            break

    return n50, l50


def calculate_gc_content(sequences):
    """
    Calculate overall GC content across all sequences

    Args:
        sequences: List of SeqRecord objects

    Returns:
        float: GC content as percentage
    """
    total_gc = 0
    total_bases = 0

    for seq in sequences:
        seq_str = str(seq.seq).upper()
        gc_count = seq_str.count('G') + seq_str.count('C')
        total_gc += gc_count
        total_bases += len(seq_str)

    if total_bases == 0:
        return 0.0

    return (total_gc / total_bases) * 100


def calculate_assembly_stats(fasta_file, sample_name=None):
    """
    Calculate comprehensive assembly statistics for a single assembly

    Args:
        fasta_file: Path to FASTA file
        sample_name: Optional sample name (for individual assemblies)

    Returns:
        dict: Assembly statistics
    """
    # Parse sequences
    sequences = parse_fasta(fasta_file)

    # Basic counts
    num_contigs = len(sequences)

    # Handle empty assemblies
    if num_contigs == 0:
        return {
            'sample': sample_name or 'coassembly',
            'num_contigs': 0,
            'total_size_bp': 0,
            'mean_length_bp': 0,
            'longest_contig_bp': 0,
            'n50': 0,
            'l50': 0,
            'gc_percent': 0.0
        }

    # Get all contig lengths
    lengths = [len(seq.seq) for seq in sequences]

    # Calculate statistics
    total_size = sum(lengths)
    mean_length = np.mean(lengths)
    longest_contig = max(lengths)
    n50, l50 = calculate_n50_l50(lengths)
    gc_content = calculate_gc_content(sequences)

    return {
        'sample': sample_name or 'coassembly',
        'num_contigs': num_contigs,
        'total_size_bp': total_size,
        'mean_length_bp': int(mean_length),
        'longest_contig_bp': longest_contig,
        'n50': n50,
        'l50': l50,
        'gc_percent': round(gc_content, 2)
    }


def main():
    """Main function called by Snakemake"""
    # Get inputs from snakemake
    contigs_input = snakemake.input.contigs
    output_file = snakemake.output.stats

    # Determine if coassembly or individual
    assembly_strategy = snakemake.config.get("pipeline", {}).get("assembly_strategy", "coassembly")

    results = []

    if assembly_strategy == "coassembly":
        # Single input file - coassembly
        if isinstance(contigs_input, str):
            stats = calculate_assembly_stats(contigs_input, sample_name="coassembly")
            results.append(stats)
        else:
            # Shouldn't happen, but handle just in case
            stats = calculate_assembly_stats(contigs_input[0], sample_name="coassembly")
            results.append(stats)

    else:
        # Multiple input files - individual assemblies
        if isinstance(contigs_input, str):
            # Single file - shouldn't happen for individual strategy
            contigs_input = [contigs_input]

        for contig_file in contigs_input:
            # Extract sample name from path
            # Expected format: .../assembly/per_sample/{sample}/final.contigs.fa
            sample_name = Path(contig_file).parent.name
            stats = calculate_assembly_stats(contig_file, sample_name=sample_name)
            results.append(stats)

    # Create DataFrame and write to TSV
    df = pd.DataFrame(results)

    # Sort by sample name
    df = df.sort_values('sample')

    # Write output
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Assembly statistics written to {output_file}")
    print(f"Processed {len(results)} assembly/assemblies")


if __name__ == '__main__':
    main()
