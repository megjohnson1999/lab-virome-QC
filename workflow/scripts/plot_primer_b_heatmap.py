"""
Visualize Primer B Cross-Contamination

Creates a heatmap showing primer B distribution across samples to identify
cross-contamination patterns.

The heatmap shows:
- Rows: Samples
- Columns: Primer B variants (3GB-1 through 3GB-24)
- Color intensity: % of primer B reads for that variant
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Snakemake inputs
distribution_file = snakemake.input.distribution
summary_file = snakemake.input.summary

# Snakemake outputs
heatmap_output = snakemake.output.heatmap

# Load data
print("Loading primer B distribution data...")
dist_df = pd.read_csv(distribution_file, sep='\t')
summary_df = pd.read_csv(summary_file, sep='\t')

# Create pivot table for heatmap (samples × primer B variants)
print("Creating heatmap matrix...")
heatmap_data = dist_df.pivot(
    index='sample',
    columns='primer_b_variant',
    values='percent_of_pb_reads'
).fillna(0)

# Sort columns naturally (3GB-1, 3GB-2, ..., 3GB-24)
def natural_sort_key(s):
    """Sort primer B names naturally (3GB-1, 3GB-2, ..., 3GB-10, 3GB-11, ...)"""
    import re
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))]

heatmap_data = heatmap_data[sorted(heatmap_data.columns, key=natural_sort_key)]

# Sort rows by dominant primer B for better visualization
sample_order = summary_df.sort_values(['dominant_pb', 'sample'])['sample'].tolist()
heatmap_data = heatmap_data.reindex(sample_order)

# Create figure with smart scaling for large datasets
n_samples = len(heatmap_data)
n_variants = len(heatmap_data.columns)

# For large datasets, focus on contaminated samples + representative set
if n_samples > 100:
    print(f"Large dataset detected ({n_samples} samples). Focusing on contaminated samples and representatives.")

    # Get contaminated samples
    contaminated_samples = summary_df[summary_df['contamination_flag']]['sample'].tolist()

    # Get a representative sample from each dominant primer B group (max 50 total)
    dominant_groups = summary_df.groupby('dominant_pb')['sample'].apply(list)
    representative_samples = []
    for pb_variant, group_samples in dominant_groups.items():
        # Take first sample from each group as representative
        if group_samples:
            representative_samples.append(group_samples[0])

    # Combine contaminated + representatives, limit to reasonable size
    focus_samples = list(set(contaminated_samples + representative_samples))
    if len(focus_samples) > 80:  # Still too many, prioritize contaminated
        focus_samples = contaminated_samples + representative_samples[:80-len(contaminated_samples)]

    # Filter heatmap data to focus samples
    heatmap_data = heatmap_data.loc[heatmap_data.index.intersection(focus_samples)]

    # Add note to title
    title_suffix = f"\nShowing {len(contaminated_samples)} contaminated + {len(focus_samples)-len(contaminated_samples)} representative samples (of {n_samples} total)"
else:
    title_suffix = ""

# Calculate reasonable figure size with limits
fig_height = max(8, min(len(heatmap_data) * 0.4, 40))  # Cap at 40 inches height
fig_width = max(12, min(n_variants * 0.5, 20))  # Cap at 20 inches width

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Create heatmap
sns.heatmap(
    heatmap_data,
    cmap='YlOrRd',
    cbar_kws={'label': '% of Primer B Reads'},
    linewidths=0.5,
    linecolor='lightgray',
    ax=ax,
    vmin=0,
    vmax=100,
    square=False
)

# Formatting
ax.set_xlabel('Primer B Variant', fontsize=12, fontweight='bold')
ax.set_ylabel('Sample', fontsize=12, fontweight='bold')
ax.set_title(f'Primer B Cross-Contamination Heatmap\n(% of Primer B Reads per Variant){title_suffix}',
             fontsize=14, fontweight='bold', pad=20)

# Rotate x-axis labels for readability
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0, fontsize=8)

# Add contamination flags as colored markers for samples shown in heatmap
contaminated_samples = summary_df[summary_df['contamination_flag']]['sample'].tolist()
for i, sample in enumerate(heatmap_data.index):
    if sample in contaminated_samples:
        # Add red marker on y-axis for contaminated samples - positioned to avoid text overlap
        ax.text(-1.2, i + 0.5, '!', ha='center', va='center',
                fontsize=14, color='red', fontweight='bold',
                bbox=dict(boxstyle="circle,pad=0.1", facecolor='red', alpha=0.8))

# Adjust layout
plt.tight_layout()

# Save figure
print(f"Saving heatmap to {heatmap_output}")
plt.savefig(heatmap_output, dpi=300, bbox_inches='tight')
plt.close()

print("\n=== Heatmap Statistics ===")
print(f"Total samples analyzed: {n_samples}")
print(f"Samples shown in heatmap: {len(heatmap_data)}")
print(f"Primer B variants: {len(heatmap_data.columns)}")
total_contaminated = len(contaminated_samples)
shown_contaminated = len([s for s in contaminated_samples if s in heatmap_data.index])
print(f"Contaminated samples (marked with red !): {shown_contaminated}/{total_contaminated} shown")

# Print dominant primer B distribution
print("\nDominant Primer B Distribution:")
dominant_counts = summary_df['dominant_pb'].value_counts()
for pb, count in dominant_counts.items():
    print(f"  {pb}: {count} samples")
