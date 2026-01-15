#!/usr/bin/env python3
"""
Generate QC quality metrics for virome samples

Evaluates samples based on 3 essential quality metrics:
- % host reads (VLP prep efficiency)
- % rRNA reads (biological contamination)
- Final read count after QC (sequencing depth)

Provides quality categories (Excellent/Good/Concerning/Poor) based on practical
thresholds rather than binary pass/fail assessments, allowing researchers to make
informed decisions about sample usability for their specific analysis goals.

"""

import pandas as pd
import sys
from pathlib import Path

# Quality assessment relies on direct contamination measurements

def categorize_quality(host_pct: float, rrna_pct: float, final_reads: int,
                      max_host_pct: float, max_rrna_pct: float, min_final_reads: int) -> str:
    """
    Categorize sample quality based on practical thresholds

    Categories:
    - Excellent: All metrics well within thresholds
    - Good: All metrics pass but some approaching thresholds
    - Concerning: One metric exceeds threshold but sample still usable
    - Poor: Multiple metrics exceed thresholds, likely problematic
    """
    host_excellent = host_pct <= max_host_pct * 0.5
    host_good = host_pct <= max_host_pct

    rrna_excellent = rrna_pct <= max_rrna_pct * 0.5
    rrna_good = rrna_pct <= max_rrna_pct

    reads_excellent = final_reads >= min_final_reads * 2
    reads_good = final_reads >= min_final_reads

    # Count issues
    issues = 0
    if not host_good:
        issues += 1
    if not rrna_good:
        issues += 1
    if not reads_good:
        issues += 1

    if issues == 0:
        if host_excellent and rrna_excellent and reads_excellent:
            return "Excellent"
        else:
            return "Good"
    elif issues == 1:
        return "Concerning"
    else:
        return "Poor"

def main():
    # Get inputs from snakemake
    read_counts_file = snakemake.input.read_counts
    output_file = snakemake.output[0]

    # Load simplified thresholds from config (3 metrics only)
    config = snakemake.config
    thresholds = config.get('qc_thresholds', {})
    max_host_pct = thresholds.get('max_host_percent', 10)
    max_rrna_pct = thresholds.get('max_rrna_percent', 20)
    min_final = thresholds.get('min_final_reads', 100000)

    # Load read counts
    read_df = pd.read_csv(read_counts_file, sep='\t')

    # Initialize results
    results = []

    # Process each sample - extract from read counts file
    samples = read_df['sample'].unique()

    for sample in samples:
        # Quality metrics structure - 3 essential metrics
        metrics = {
            'sample': sample,
            'host_percent': 'NA',
            'rrna_percent': 'NA',
            'final_reads': 'NA',
            'quality_category': 'Unknown',
            'notes': []
        }

        # Calculate % reads retained at each step
        sample_reads = read_df[read_df['sample'] == sample]

        if len(sample_reads) > 0:
            fastp_count = sample_reads[sample_reads['step'] == 'fastp']['reads'].values[0]
            host_depleted = sample_reads[sample_reads['step'] == 'host_depleted']['reads'].values[0]
            clean_count = sample_reads[sample_reads['step'] == 'clean']['reads'].values[0]

            # Calculate % host removed (using fastp count as denominator)
            # fastp is the step immediately before host depletion, so this gives
            # the true host contamination % without including QC filtering losses
            host_removed = fastp_count - host_depleted
            host_pct = (host_removed / fastp_count) * 100 if fastp_count > 0 else 0
            metrics['host_percent'] = host_pct

            # Calculate % pure rRNA contamination (biological contamination only)
            # Use the actual rrna_removed step count rather than combined processing losses
            rrna_removed_step = sample_reads[sample_reads['step'] == 'rrna_removed']['reads']
            if len(rrna_removed_step) > 0:
                rrna_removed_count = rrna_removed_step.values[0]
                # Pure rRNA contamination: reads lost during rRNA removal step only
                rrna_contamination = host_depleted - rrna_removed_count
                rrna_pct = (rrna_contamination / host_depleted) * 100 if host_depleted > 0 else 0
            else:
                # Fallback to old method if rrna_removed step not found
                rrna_contamination = host_depleted - clean_count
                rrna_pct = (rrna_contamination / host_depleted) * 100 if host_depleted > 0 else 0
            metrics['rrna_percent'] = rrna_pct
            metrics['final_reads'] = clean_count

            # Categorize overall quality
            metrics['quality_category'] = categorize_quality(
                host_pct, rrna_pct, clean_count,
                max_host_pct, max_rrna_pct, min_final
            )

            # Add informational notes (not failure indicators)
            if host_pct > max_host_pct:
                metrics['notes'].append(f'High_host({host_pct:.1f}%)')
            if rrna_pct > max_rrna_pct:
                metrics['notes'].append(f'High_rRNA({rrna_pct:.1f}%)')
            if clean_count < min_final:
                metrics['notes'].append(f'Low_coverage({int(clean_count)})')

        metrics['notes'] = ';'.join(metrics['notes']) if metrics['notes'] else 'None'
        results.append(metrics)

    # Write results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, sep='\t', index=False)

    # Print summary
    print("\n" + "="*60)
    print("VIROME QC METRICS SUMMARY (3 Essential Metrics)")
    print("="*60)
    print(f"Samples evaluated: {len(results_df)}")

    # Count samples by quality category
    quality_counts = results_df['quality_category'].value_counts()
    for category in ['Excellent', 'Good', 'Concerning', 'Poor', 'Unknown']:
        count = quality_counts.get(category, 0)
        print(f"{category}: {count}")

    print("")
    print("Sample Details:")
    print("-" * 60)
    print(results_df.to_string(index=False))
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
