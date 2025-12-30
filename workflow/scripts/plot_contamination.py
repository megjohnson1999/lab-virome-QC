#!/usr/bin/env python3
"""
Generate contamination visualization plots

Creates publication-quality plots emphasizing outlier detection:
1. Bar plot highlighting samples that deviate from the group
2. Box plots showing distribution with outliers marked
3. Scatter plot showing correlation with outliers labeled
4. Heatmap showing patterns across samples
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path


def identify_outliers(series):
    """
    Identify statistical outliers using IQR method
    Returns boolean mask of outliers
    """
    if len(series) < 4:
        return pd.Series([False] * len(series), index=series.index)

    Q1 = series.quantile(0.25)
    Q3 = series.quantile(0.75)
    IQR = Q3 - Q1

    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    return (series < lower_bound) | (series > upper_bound)


def plot_contamination_bars(df, output_prefix):
    """
    Create stacked bar plot with clean professional styling to match rRNA histogram
    """
    # Dynamic sizing based on sample count
    n_samples = len(df)
    fig_width = max(16, n_samples * 0.7)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    # Identify outliers
    df['is_outlier'] = identify_outliers(df['total_contamination_percent'])

    # Sort by total contamination
    df_sorted = df.sort_values('total_contamination_percent', ascending=True).copy()

    # Create bars with steel blue theme (matching rRNA plot)
    x = range(len(df_sorted))
    width = 0.8

    # Use steel blue color scheme for consistency
    phix_color = '#5B9BD5'  # Steel blue (matching rRNA)
    vector_color = '#8DB4E3'  # Lighter steel blue for vector
    outlier_phix_color = '#DC3545'  # Red for outliers
    outlier_vector_color = '#F8D7DA'  # Light red for outlier vectors

    # Color outliers differently
    phix_colors = [outlier_phix_color if outlier else phix_color for outlier in df_sorted['is_outlier']]
    vector_colors = [outlier_vector_color if outlier else vector_color for outlier in df_sorted['is_outlier']]

    # Plot bars with clean edges
    ax.bar(x, df_sorted['phix_percent'], width, color=phix_colors,
           edgecolor='#2F5F8F', linewidth=0.8, label='PhiX')
    ax.bar(x, df_sorted['vector_percent'], width,
           bottom=df_sorted['phix_percent'],
           color=vector_colors, edgecolor='#2F5F8F', linewidth=0.8, label='Vector/Plasmid')

    # Calculate statistics for reference lines
    mean_total = df['total_contamination_percent'].mean()
    median_total = df['total_contamination_percent'].median()

    # Add mean and median reference lines (matching rRNA style)
    ax.axhline(y=mean_total, color='#FF8C00', linestyle='--', linewidth=2.5,
               alpha=0.9, label=f'Mean: {mean_total:.3f}%')
    ax.axhline(y=median_total, color='#DC3545', linestyle='--', linewidth=2.5,
               alpha=0.9, label=f'Median: {median_total:.3f}%')

    # Formatting with Arial font (matching rRNA plot)
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.set_ylabel('Contamination (%)', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.set_title(f'Contamination Distribution (n={n_samples})\nPhiX174 + Vector/Plasmid Detection',
                fontsize=14, fontweight='bold', fontfamily='Arial')
    ax.set_xticks(x)

    # Adjust font size based on sample count
    label_fontsize = 8 if n_samples > 40 else 9 if n_samples > 30 else 10
    ax.set_xticklabels(df_sorted['sample'], rotation=60, ha='right', fontsize=label_fontsize)

    # Highlight outlier sample names in red
    for i, (idx, row) in enumerate(df_sorted.iterrows()):
        if row['is_outlier']:
            ax.get_xticklabels()[i].set_color('#DC3545')
            ax.get_xticklabels()[i].set_fontweight('bold')

    # Clean grid (matching rRNA style)
    ax.grid(axis='y', alpha=0.4, linestyle='-', linewidth=0.5, color='#E5E5E5')
    ax.set_facecolor('white')

    # Legend with clean styling
    ax.legend(loc='upper left', frameon=True, fancybox=False, fontsize=10,
              edgecolor='#E5E5E5', framealpha=0.95)

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_bars.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved bar plot: {output_file}")


def plot_contamination_boxes(df, output_prefix):
    """
    Create box plots with clean professional styling to match rRNA histogram
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 6))

    # Identify outliers for each type
    phix_outliers = identify_outliers(df['phix_percent'])
    vector_outliers = identify_outliers(df['vector_percent'])
    total_outliers = identify_outliers(df['total_contamination_percent'])

    # Steel blue color scheme (matching rRNA plot)
    phix_color = '#5B9BD5'
    vector_color = '#8DB4E3'
    total_color = '#A5C8E7'

    # PhiX distribution
    bp1 = axes[0].boxplot([df['phix_percent']], labels=['PhiX174'],
                          patch_artist=True,
                          boxprops=dict(facecolor=phix_color, alpha=0.8, linewidth=1.5),
                          medianprops=dict(color='#2F5F8F', linewidth=2.5),
                          whiskerprops=dict(linewidth=1.5, color='#2F5F8F'),
                          capprops=dict(linewidth=1.5, color='#2F5F8F'),
                          flierprops=dict(marker='o', markerfacecolor='#DC3545', markersize=6,
                                        markeredgecolor='#DC3545', alpha=0.9))
    axes[0].set_ylabel('Contamination (%)', fontsize=12, fontweight='bold', fontfamily='Arial')
    axes[0].set_title(f'PhiX174 Distribution\n({phix_outliers.sum()} outliers)',
                      fontsize=12, fontweight='bold', fontfamily='Arial')
    axes[0].grid(axis='y', alpha=0.4, linestyle='-', linewidth=0.5, color='#E5E5E5')
    axes[0].set_facecolor('white')

    # Vector distribution
    bp2 = axes[1].boxplot([df['vector_percent']], labels=['Vector/Plasmid'],
                          patch_artist=True,
                          boxprops=dict(facecolor=vector_color, alpha=0.8, linewidth=1.5),
                          medianprops=dict(color='#2F5F8F', linewidth=2.5),
                          whiskerprops=dict(linewidth=1.5, color='#2F5F8F'),
                          capprops=dict(linewidth=1.5, color='#2F5F8F'),
                          flierprops=dict(marker='o', markerfacecolor='#DC3545', markersize=6,
                                        markeredgecolor='#DC3545', alpha=0.9))
    axes[1].set_ylabel('Contamination (%)', fontsize=12, fontweight='bold', fontfamily='Arial')
    axes[1].set_title(f'Vector/Plasmid Distribution\n({vector_outliers.sum()} outliers)',
                      fontsize=12, fontweight='bold', fontfamily='Arial')
    axes[1].grid(axis='y', alpha=0.4, linestyle='-', linewidth=0.5, color='#E5E5E5')
    axes[1].set_facecolor('white')

    # Total contamination distribution
    bp3 = axes[2].boxplot([df['total_contamination_percent']], labels=['Total'],
                          patch_artist=True,
                          boxprops=dict(facecolor=total_color, alpha=0.8, linewidth=1.5),
                          medianprops=dict(color='#2F5F8F', linewidth=2.5),
                          whiskerprops=dict(linewidth=1.5, color='#2F5F8F'),
                          capprops=dict(linewidth=1.5, color='#2F5F8F'),
                          flierprops=dict(marker='o', markerfacecolor='#DC3545', markersize=6,
                                        markeredgecolor='#DC3545', alpha=0.9))
    axes[2].set_ylabel('Contamination (%)', fontsize=12, fontweight='bold', fontfamily='Arial')
    axes[2].set_title(f'Total Contamination Distribution\n({total_outliers.sum()} outliers)',
                      fontsize=12, fontweight='bold', fontfamily='Arial')
    axes[2].grid(axis='y', alpha=0.4, linestyle='-', linewidth=0.5, color='#E5E5E5')
    axes[2].set_facecolor('white')

    # Add statistics text with clean styling
    for ax, data, name in zip(axes,
                               [df['phix_percent'], df['vector_percent'], df['total_contamination_percent']],
                               ['PhiX', 'Vector', 'Total']):
        median_val = data.median()
        mean_val = data.mean()
        iqr_val = data.quantile(0.75) - data.quantile(0.25)

        stats_text = f"Mean: {mean_val:.3f}%\nMedian: {median_val:.3f}%\nIQR: {iqr_val:.3f}%"
        ax.text(0.97, 0.97, stats_text, transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#E5E5E5'),
               fontsize=9, fontfamily='Arial')

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_boxes.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved box plot: {output_file}")


def plot_contamination_scatter(df, output_prefix):
    """
    Create scatter plot with clean professional styling to match rRNA histogram
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Identify outliers
    phix_outliers = identify_outliers(df['phix_percent'])
    vector_outliers = identify_outliers(df['vector_percent'])
    any_outlier = phix_outliers | vector_outliers

    # Separate outliers from normal samples
    df_normal = df[~any_outlier]
    df_outliers = df[any_outlier]

    # Plot normal samples with steel blue theme
    if len(df_normal) > 0:
        ax.scatter(df_normal['phix_percent'], df_normal['vector_percent'],
                  s=120, alpha=0.7, c='#5B9BD5',
                  edgecolors='#2F5F8F', linewidth=1,
                  label='Within expected range')

    # Plot outliers with red coloring
    if len(df_outliers) > 0:
        ax.scatter(df_outliers['phix_percent'], df_outliers['vector_percent'],
                   s=150, alpha=0.9, c='#DC3545',
                   edgecolors='#A02331', linewidth=2,
                   marker='D', label='Statistical outliers')

    # Calculate statistics for reference lines
    phix_mean = df['phix_percent'].mean()
    phix_median = df['phix_percent'].median()
    vector_mean = df['vector_percent'].mean()
    vector_median = df['vector_percent'].median()

    # Add mean and median reference lines (matching rRNA style)
    ax.axhline(y=vector_mean, color='#FF8C00', linestyle='--', linewidth=2,
               alpha=0.8, label=f'Vector Mean: {vector_mean:.3f}%')
    ax.axhline(y=vector_median, color='#DC3545', linestyle='--', linewidth=2,
               alpha=0.8, label=f'Vector Median: {vector_median:.3f}%')
    ax.axvline(x=phix_mean, color='#FF8C00', linestyle='--', linewidth=2,
               alpha=0.8, label=f'PhiX Mean: {phix_mean:.3f}%')
    ax.axvline(x=phix_median, color='#DC3545', linestyle='--', linewidth=2,
               alpha=0.8, label=f'PhiX Median: {phix_median:.3f}%')

    # Label outliers with simplified annotation
    if len(df_outliers) > 0:
        for _, row in df_outliers.iterrows():
            # Simple annotation without overcrowding
            sample_short = row['sample'].split('_')[-1] if '_' in row['sample'] else row['sample']
            if len(sample_short) > 10:
                sample_short = sample_short[:10] + '...'

            ax.annotate(sample_short,
                       (row['phix_percent'], row['vector_percent']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, fontweight='bold', color='#DC3545',
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                alpha=0.8, edgecolor='#DC3545'),
                       ha='left')

    # Formatting with Arial font (matching rRNA plot)
    ax.set_xlabel('PhiX174 Contamination (%)', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.set_ylabel('Vector/Plasmid Contamination (%)', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.set_title(f'Contamination Correlation (n={len(df)})\nPhiX174 vs Vector/Plasmid Detection',
                fontsize=14, fontweight='bold', fontfamily='Arial')

    # Clean grid (matching rRNA style)
    ax.grid(True, alpha=0.4, linestyle='-', linewidth=0.5, color='#E5E5E5')
    ax.set_facecolor('white')

    # Legend with clean styling
    ax.legend(loc='upper right', frameon=True, fancybox=False, fontsize=9,
              edgecolor='#E5E5E5', framealpha=0.95)

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_scatter.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved scatter plot: {output_file}")


def plot_contamination_heatmap(df, output_prefix):
    """
    Create heatmap with clean professional styling to match rRNA histogram
    """
    # Identify outliers
    total_outliers = identify_outliers(df['total_contamination_percent'])

    # Prepare data for heatmap
    heatmap_data = df[['sample', 'phix_percent', 'vector_percent']].set_index('sample').T

    # Sort columns by total contamination
    col_order = df.sort_values('total_contamination_percent', ascending=False)['sample']
    heatmap_data = heatmap_data[col_order]

    # Dynamic sizing based on sample count
    n_samples = len(df)
    fig_width = max(16, n_samples * 0.5)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    # Create heatmap with steel blue color scheme (similar to rRNA styling)
    # Using a colormap that matches the steel blue theme
    sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='Blues',
                cbar_kws={'label': 'Contamination (%)', 'shrink': 0.8},
                linewidths=0.8, linecolor='white', ax=ax,
                vmin=0, vmax=df[['phix_percent', 'vector_percent']].max().max(),
                square=False)

    # Highlight outlier columns with red outline (matching rRNA outlier style)
    outlier_samples = df[total_outliers]['sample'].values
    for i, sample in enumerate(col_order):
        if sample in outlier_samples:
            # Add red box around outlier columns
            col_idx = i
            ax.add_patch(plt.Rectangle((col_idx, 0), 1, 2, fill=False,
                                      edgecolor='#DC3545', lw=3, zorder=10))

    # Formatting with Arial font (matching rRNA plot)
    ax.set_xlabel('Sample (sorted by total contamination)', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.set_ylabel('Contamination Type', fontsize=12, fontweight='bold', fontfamily='Arial')
    ax.set_yticklabels(['PhiX174', 'Vector/Plasmid'], rotation=0, fontsize=11, fontfamily='Arial')

    # Adjust x-axis label rotation and font size based on sample count
    label_fontsize = 7 if n_samples > 40 else 8 if n_samples > 30 else 9
    xticklabels = ax.get_xticklabels()
    ax.set_xticklabels(xticklabels, rotation=90, fontsize=label_fontsize, fontfamily='Arial')

    # Color outlier sample names in red
    xticklabels = ax.get_xticklabels()
    for i, label in enumerate(xticklabels):
        if label.get_text() in outlier_samples:
            label.set_color('#DC3545')
            label.set_fontweight('bold')

    # Title with sample count (matching rRNA style)
    ax.set_title(f'Contamination Heatmap (n={n_samples})\nPhiX174 + Vector/Plasmid Detection',
                fontsize=14, fontweight='bold', fontfamily='Arial')

    # Add statistics annotation with clean styling
    median_total = df['total_contamination_percent'].median()
    mean_total = df['total_contamination_percent'].mean()
    outlier_count = total_outliers.sum()

    stats_text = f"Mean: {mean_total:.3f}%\nMedian: {median_total:.3f}%\nOutliers: {outlier_count}"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
           verticalalignment='top', horizontalalignment='left',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#E5E5E5'),
           fontsize=10, fontweight='bold', fontfamily='Arial')

    plt.tight_layout()

    # Save
    output_file = f"{output_prefix}_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved heatmap: {output_file}")


def main():
    # Get inputs from snakemake
    contamination_table = snakemake.input[0]
    output_prefix = snakemake.params.output_prefix

    # Load data
    df = pd.read_csv(contamination_table, sep='\t')

    # Set clean style to match rRNA histogram
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = ['Arial', 'sans-serif']
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.edgecolor'] = '#E5E5E5'
    plt.rcParams['grid.color'] = '#E5E5E5'
    plt.rcParams['grid.alpha'] = 0.4

    # Generate plots
    print("\n" + "="*60)
    print("GENERATING CONTAMINATION PLOTS WITH PROFESSIONAL STYLING")
    print("="*60 + "\n")

    plot_contamination_bars(df, output_prefix)
    plot_contamination_boxes(df, output_prefix)
    plot_contamination_scatter(df, output_prefix)
    plot_contamination_heatmap(df, output_prefix)

    print("\n" + "="*60)
    print("CONTAMINATION PLOTS COMPLETE - STYLED TO MATCH rRNA HISTOGRAM")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
