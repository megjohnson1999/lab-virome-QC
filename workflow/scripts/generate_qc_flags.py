#!/usr/bin/env python3
"""
Generate QC pass/fail flags for virome samples

Evaluates samples based on 3 essential quality metrics:
- % host reads (VLP prep efficiency)
- % rRNA reads (biological contamination)
- Final read count after QC (sequencing depth)

Note: ViromeQC enrichment scoring removed as redundant with direct contamination measurements
"""

import pandas as pd
import sys
from pathlib import Path

# ViromeQC parsing function removed - enrichment scoring no longer used
# Quality assessment now relies on direct contamination measurements

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
        # Simplified flags structure - 3 quality metrics only
        flags = {
            'sample': sample,
            'host_percent': 'NA',
            'pass_host': 'UNKNOWN',
            'rrna_percent': 'NA',
            'pass_rrna': 'UNKNOWN',
            'pass_final_count': 'UNKNOWN',
            'overall_pass': 'UNKNOWN',
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
            flags['host_percent'] = host_pct

            if host_pct <= max_host_pct:
                flags['pass_host'] = 'PASS'
            else:
                flags['pass_host'] = 'FAIL'
                flags['notes'].append(f'High_host({host_pct:.1f}%)')

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
            flags['rrna_percent'] = rrna_pct

            if rrna_pct <= max_rrna_pct:
                flags['pass_rrna'] = 'PASS'
            else:
                flags['pass_rrna'] = 'FAIL'
                flags['notes'].append(f'High_rRNA({rrna_pct:.1f}%)')

            # Check final read count
            if clean_count >= min_final:
                flags['pass_final_count'] = 'PASS'
            else:
                flags['pass_final_count'] = 'FAIL'
                flags['notes'].append(f'Low_final_reads({int(clean_count)})')

        # Overall pass/fail - simplified to 3 essential metrics
        all_checks = [
            flags['pass_host'],
            flags['pass_rrna'],
            flags['pass_final_count']
        ]

        if 'FAIL' in all_checks:
            flags['overall_pass'] = 'FAIL'
        elif 'UNKNOWN' in all_checks:
            flags['overall_pass'] = 'WARNING'
        else:
            flags['overall_pass'] = 'PASS'

        flags['notes'] = ';'.join(flags['notes']) if flags['notes'] else 'None'
        results.append(flags)

    # Write results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, sep='\t', index=False)

    # Print summary
    print("\n" + "="*60)
    print("SIMPLIFIED QC FLAGS SUMMARY (3 Essential Metrics)")
    print("="*60)
    print(f"Samples evaluated: {len(results_df)}")
    print(f"PASS: {len(results_df[results_df['overall_pass'] == 'PASS'])}")
    print(f"FAIL: {len(results_df[results_df['overall_pass'] == 'FAIL'])}")
    print(f"WARNING: {len(results_df[results_df['overall_pass'] == 'WARNING'])}")
    print("")
    print(results_df.to_string(index=False))
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
