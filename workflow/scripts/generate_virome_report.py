#!/usr/bin/env python3
"""
Generate user-friendly virome QC report (HTML + PDF)

This script consolidates all QC data sources into a unified structure,
implements enhanced outlier detection, and generates both interactive HTML
dashboard and publication-ready PDF reports.

Replaces MultiQC with virome-specific reporting focused on:
- Batch-level overview and outlier identification
- Clear pass/fail assessment
- Sample ranking by quality metrics
- Interactive exploration and static sharing
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
import re
import sys
from typing import Dict, List, Tuple, Optional, Any
import warnings
from datetime import datetime

# HTML template generation
try:
    from jinja2 import Environment, FileSystemLoader, Template
    JINJA2_AVAILABLE = True
except ImportError:
    JINJA2_AVAILABLE = False
    print("Warning: jinja2 not available. HTML generation will be skipped.")

# Suppress pandas warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)

class ViromeQCDataLoader:
    """
    Unified data loader for all virome QC data sources
    """

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.qc_thresholds = config.get('qc_thresholds', {})
        self.samples = []
        self.data = {}

    def load_read_counts(self, read_counts_file: str) -> pd.DataFrame:
        """Load and process read count data"""
        read_df = pd.read_csv(read_counts_file, sep='\t')

        # Pivot to have samples as rows and steps as columns
        read_pivot = read_df.pivot(index='sample', columns='step', values='reads')

        # Calculate retention percentages
        if 'raw' in read_pivot.columns:
            for col in read_pivot.columns:
                if col != 'raw':
                    read_pivot[f'{col}_retention'] = (read_pivot[col] / read_pivot['raw'] * 100)

        return read_pivot

    def load_contamination_data(self, contamination_file: str) -> pd.DataFrame:
        """Load contamination summary data"""
        return pd.read_csv(contamination_file, sep='\t')

    def load_qc_flags(self, qc_flags_file: str) -> pd.DataFrame:
        """Load QC pass/fail flags"""
        return pd.read_csv(qc_flags_file, sep='\t')

    def load_primer_b_data(self, primer_b_file: str) -> pd.DataFrame:
        """Load primer B contamination data"""
        return pd.read_csv(primer_b_file, sep='\t')

    def load_viromeqc_data(self, viromeqc_files: List[str]) -> pd.DataFrame:
        """Load and parse ViromeQC enrichment scores"""
        viromeqc_data = []

        for vqc_file in viromeqc_files:
            # Extract sample name from file path
            sample = Path(vqc_file).stem.replace('_viromeqc', '')

            enrichment_score = None
            with open(vqc_file) as f:
                for line in f:
                    if "Enrichment score" in line or "enrichment" in line.lower():
                        parts = line.strip().split()
                        try:
                            enrichment_score = float(parts[-1])
                            break
                        except (ValueError, IndexError):
                            pass

            viromeqc_data.append({
                'sample': sample,
                'enrichment_score': enrichment_score
            })

        return pd.DataFrame(viromeqc_data)

    def parse_fastqc_summary(self, fastqc_files: List[str]) -> pd.DataFrame:
        """Parse FastQC summary statistics from ZIP files"""
        fastqc_data = []

        for fastqc_file in fastqc_files:
            # Extract sample name and stage (raw/trimmed/final)
            file_path = Path(fastqc_file)
            filename = file_path.stem

            # Parse sample name and stage from filename
            # Format: {sample}_R1_fastqc.zip in raw/trimmed/final directories
            stage = file_path.parent.name  # raw, trimmed, or final
            sample = filename.replace('_R1_fastqc', '').replace('_R2_fastqc', '')

            # For now, store basic metadata - can extend to parse ZIP contents if needed
            fastqc_data.append({
                'sample': sample,
                'stage': stage,
                'file_path': str(fastqc_file)
            })

        return pd.DataFrame(fastqc_data)

class ViromeQCOutlierDetector:
    """
    Enhanced outlier detection across all QC metrics
    """

    def __init__(self):
        self.outlier_methods = {
            'iqr': self._iqr_outliers,
            'zscore': self._zscore_outliers,
            'modified_zscore': self._modified_zscore_outliers
        }

    def _iqr_outliers(self, series: pd.Series) -> pd.Series:
        """IQR-based outlier detection (existing method)"""
        if len(series) < 4:
            return pd.Series([False] * len(series), index=series.index)

        Q1 = series.quantile(0.25)
        Q3 = series.quantile(0.75)
        IQR = Q3 - Q1

        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        return (series < lower_bound) | (series > upper_bound)

    def _zscore_outliers(self, series: pd.Series, threshold: float = 2.0) -> pd.Series:
        """Z-score based outlier detection"""
        if len(series) < 4:
            return pd.Series([False] * len(series), index=series.index)

        z_scores = np.abs((series - series.mean()) / series.std())
        return z_scores > threshold

    def _modified_zscore_outliers(self, series: pd.Series, threshold: float = 3.5) -> pd.Series:
        """Modified Z-score using median absolute deviation"""
        if len(series) < 4:
            return pd.Series([False] * len(series), index=series.index)

        median = series.median()
        mad = np.median(np.abs(series - median))

        if mad == 0:
            return pd.Series([False] * len(series), index=series.index)

        modified_z_scores = 0.6745 * (series - median) / mad
        return np.abs(modified_z_scores) > threshold

    def detect_outliers(self, data: pd.DataFrame, columns: List[str],
                       method: str = 'iqr') -> pd.DataFrame:
        """
        Detect outliers across multiple metrics
        """
        outlier_df = data.copy()

        if method not in self.outlier_methods:
            raise ValueError(f"Unknown outlier method: {method}")

        detect_func = self.outlier_methods[method]

        for col in columns:
            if col in data.columns:
                outlier_df[f'{col}_outlier'] = detect_func(data[col])

        return outlier_df

class ViromeQCReportGenerator:
    """
    Main class for generating virome QC reports
    """

    def __init__(self, inputs: Dict[str, Any], outputs: Dict[str, str],
                 config: Dict[str, Any]):
        self.inputs = inputs
        self.outputs = outputs
        self.config = config

        self.loader = ViromeQCDataLoader(config)
        self.outlier_detector = ViromeQCOutlierDetector()
        self.unified_data = None

    def load_all_data(self):
        """Load and unify all QC data sources"""
        print("Loading QC data sources...")

        # Load individual data sources
        read_counts = self.loader.load_read_counts(self.inputs['read_counts'])
        contamination = self.loader.load_contamination_data(self.inputs['contamination_summary'])
        qc_flags = self.loader.load_qc_flags(self.inputs['qc_flags'])
        primer_b = self.loader.load_primer_b_data(self.inputs['primer_b_summary'])
        viromeqc = self.loader.load_viromeqc_data(self.inputs['viromeqc_files'])

        # Start with QC flags as the base (has all samples)
        unified = qc_flags.set_index('sample')

        # Merge other data sources
        if not read_counts.empty:
            unified = unified.join(read_counts, how='left', rsuffix='_reads')

        if not contamination.empty:
            contamination_indexed = contamination.set_index('sample')
            unified = unified.join(contamination_indexed, how='left', rsuffix='_contam')

        if not primer_b.empty:
            primer_b_indexed = primer_b.set_index('sample')
            unified = unified.join(primer_b_indexed, how='left', rsuffix='_pb')

        if not viromeqc.empty:
            viromeqc_indexed = viromeqc.set_index('sample')
            unified = unified.join(viromeqc_indexed, how='left', rsuffix='_vqc')

        # Reset index to have sample as column
        unified = unified.reset_index()

        self.unified_data = unified
        print(f"Successfully loaded data for {len(unified)} samples")

    def calculate_quality_scores(self):
        """Calculate composite quality scores for sample ranking"""
        if self.unified_data is None:
            raise ValueError("Must load data first")

        # Define quality metrics and their weights
        quality_metrics = {
            'enrichment_score': 0.3,  # ViromeQC enrichment (higher = better)
            'clean_retention': 0.25,   # Final read retention (higher = better)
            'phix_percent': -0.2,      # PhiX contamination (lower = better)
            'vector_percent': -0.15,   # Vector contamination (lower = better)
            'total_contamination_percent': -0.1  # Total contamination (lower = better)
        }

        # Normalize metrics to 0-1 scale
        quality_scores = []

        for idx, row in self.unified_data.iterrows():
            score = 0.0
            total_weight = 0.0

            for metric, weight in quality_metrics.items():
                if metric in row and pd.notna(row[metric]):
                    value = row[metric]

                    if metric == 'enrichment_score':
                        # Normalize enrichment score (cap at 50 for normalization)
                        normalized = min(value, 50) / 50
                    elif metric.endswith('_retention'):
                        # Retention percentage (already 0-100)
                        normalized = value / 100
                    elif metric.endswith('_percent'):
                        # Contamination percentage (invert since lower is better)
                        # Cap at 10% for normalization
                        normalized = max(0, 1 - min(value, 10) / 10)
                    else:
                        normalized = value

                    score += weight * normalized
                    total_weight += abs(weight)

            # Normalize by total weight used
            if total_weight > 0:
                quality_scores.append(score / total_weight * 100)  # Convert to 0-100 scale
            else:
                quality_scores.append(0)

        self.unified_data['quality_score'] = quality_scores

        # Rank samples by quality score
        self.unified_data['quality_rank'] = self.unified_data['quality_score'].rank(
            method='dense', ascending=False
        ).astype(int)

    def detect_all_outliers(self):
        """Detect outliers across all quality metrics"""
        outlier_columns = [
            'enrichment_score',
            'clean_retention',
            'phix_percent',
            'vector_percent',
            'total_contamination_percent',
            'quality_score'
        ]

        # Filter to columns that exist in the data
        existing_columns = [col for col in outlier_columns if col in self.unified_data.columns]

        self.unified_data = self.outlier_detector.detect_outliers(
            self.unified_data, existing_columns, method='iqr'
        )

        # Add composite outlier flag (outlier in any metric)
        outlier_cols = [col for col in self.unified_data.columns if col.endswith('_outlier')]
        self.unified_data['any_outlier'] = self.unified_data[outlier_cols].any(axis=1)

    def calculate_batch_statistics(self) -> Dict[str, Any]:
        """Calculate batch-level summary statistics"""
        stats = {
            'total_samples': len(self.unified_data),
            'pass_count': len(self.unified_data[self.unified_data['overall_pass'] == 'PASS']),
            'fail_count': len(self.unified_data[self.unified_data['overall_pass'] == 'FAIL']),
            'warning_count': len(self.unified_data[self.unified_data['overall_pass'] == 'WARNING']),
            'outlier_count': len(self.unified_data[self.unified_data['any_outlier'] == True]),
        }

        stats['pass_rate'] = stats['pass_count'] / stats['total_samples'] * 100

        # Calculate median values for key metrics
        numeric_cols = ['enrichment_score', 'clean_retention', 'phix_percent',
                       'vector_percent', 'quality_score']

        for col in numeric_cols:
            if col in self.unified_data.columns:
                stats[f'{col}_median'] = self.unified_data[col].median()
                stats[f'{col}_mean'] = self.unified_data[col].mean()
                stats[f'{col}_std'] = self.unified_data[col].std()

        return stats

    def generate_json_data(self) -> Dict[str, Any]:
        """Generate structured data for HTML/PDF rendering"""
        batch_stats = self.calculate_batch_statistics()

        # Convert DataFrame to dict for JSON serialization
        samples_data = self.unified_data.to_dict('records')

        # Get top outliers for highlighting
        outlier_samples = self.unified_data[self.unified_data['any_outlier'] == True]
        top_outliers = outlier_samples.nsmallest(5, 'quality_score').to_dict('records')

        return {
            'batch_statistics': batch_stats,
            'samples': samples_data,
            'top_outliers': top_outliers,
            'config': self.config,
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }

    def generate_html_report(self, template_dir: str, output_file: str) -> bool:
        """
        Generate HTML report using Jinja2 template
        """
        if not JINJA2_AVAILABLE:
            print("Skipping HTML generation - jinja2 not available")
            return False

        try:
            # Set up Jinja2 environment
            env = Environment(loader=FileSystemLoader(template_dir))
            template = env.get_template('virome_report.html')

            # Generate report data
            report_data = self.generate_json_data()

            # Render template
            html_content = template.render(
                batch_statistics=report_data['batch_statistics'],
                samples=report_data['samples'],
                top_outliers=report_data['top_outliers'],
                config=report_data['config'],
                timestamp=report_data['timestamp'],
                report_data=report_data  # Full data for JavaScript
            )

            # Write HTML file
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html_content)

            print(f"HTML report generated: {output_file}")
            return True

        except Exception as e:
            print(f"Error generating HTML report: {e}")
            import traceback
            traceback.print_exc()
            return False

    def generate(self):
        """Main generation method"""
        print("="*60)
        print("VIROME QC REPORT GENERATOR")
        print("="*60)

        # Phase 1: Load and unify all data
        self.load_all_data()

        # Phase 2: Calculate quality metrics
        self.calculate_quality_scores()

        # Phase 3: Detect outliers
        self.detect_all_outliers()

        # Phase 4: Generate output data structure
        report_data = self.generate_json_data()

        # Phase 5: Save structured data
        output_json = self.outputs.get('json_data', 'virome_report_data.json')
        with open(output_json, 'w') as f:
            json.dump(report_data, f, indent=2, default=str)

        print(f"\nGenerated structured report data: {output_json}")

        # Phase 6: Generate HTML report
        html_output = self.outputs.get('html_report')
        if html_output:
            # Determine template directory relative to script location
            script_dir = Path(__file__).parent
            template_dir = script_dir.parent / 'templates'

            success = self.generate_html_report(str(template_dir), html_output)
            if success:
                print(f"HTML report generated: {html_output}")
            else:
                print("Failed to generate HTML report")

        # Phase 7: Generate PDF report (placeholder for future implementation)
        pdf_output = self.outputs.get('pdf_report')
        if pdf_output:
            print(f"PDF generation not yet implemented, skipping: {pdf_output}")

        # Print summary
        stats = report_data['batch_statistics']
        print(f"\nBATCH SUMMARY:")
        print(f"  Total samples: {stats['total_samples']}")
        print(f"  Pass rate: {stats['pass_rate']:.1f}%")
        print(f"  Outliers detected: {stats['outlier_count']}")
        print(f"  Quality score median: {stats.get('quality_score_median', 'N/A'):.1f}")

        print("\n" + "="*60)
        return report_data

def main():
    """Main entry point for Snakemake integration"""
    # Get inputs and outputs from snakemake
    inputs = {
        'read_counts': snakemake.input.read_counts,
        'contamination_summary': snakemake.input.contamination_summary,
        'qc_flags': snakemake.input.qc_flags,
        'primer_b_summary': snakemake.input.primer_b_summary,
        'viromeqc_files': snakemake.input.viromeqc_files,
        'fastqc_files': snakemake.input.get('fastqc_files', [])
    }

    outputs = {
        'json_data': snakemake.output.json_data,
        'html_report': snakemake.output.get('html_report'),
        'pdf_report': snakemake.output.get('pdf_report')
    }

    config = snakemake.config

    # Generate report
    generator = ViromeQCReportGenerator(inputs, outputs, config)
    report_data = generator.generate()

    return report_data

if __name__ == "__main__":
    main()