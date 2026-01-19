#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 10:40:41 2026

@author: alicia-pageau
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import math
from matplotlib import colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FixedLocator
from matplotlib.font_manager import FontProperties

wkdir = '/home/alicia-pageau/Documents/antifungal_project/PDR1/00_scripts/Jann_et_al_2025/Upadated_scripts_after_review_11_2025/Figure1/'
os.chdir(wkdir)

# aa order per properties for heatmaps
aa_order = ["*", "P", "G", "C", "Q", "N", "T", "S", "E", "D",
            "K", "H", "R", "W", "Y", "F", "M", "L", "I", "V", "A"]

codon_to_aa = {
    'TTT': 'F', 'TTC': 'F',  # Phenylalanine
    'TTA': 'L', 'TTG': 'L',  # Leucine
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',  # Leucine
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',  # Isoleucine
    'ATG': 'M',  # Methionine (start codon)
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',  # Valine
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',  # Serine
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',  # Proline
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',  # Threonine
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',  # Alanine
    'TAT': 'Y', 'TAC': 'Y',  # Tyrosine
    'TAA': '*', 'TAG': '*',  # Stop codons
    'CAT': 'H', 'CAC': 'H',  # Histidine
    'CAA': 'Q', 'CAG': 'Q',  # Glutamine
    'AAT': 'N', 'AAC': 'N',  # Asparagine
    'AAA': 'K', 'AAG': 'K',  # Lysine
    'GAT': 'D', 'GAC': 'D',  # Aspartic acid
    'GAA': 'E', 'GAG': 'E',  # Glutamic acid
    'TGT': 'C', 'TGC': 'C',  # Cysteine
    'TGA': '*',  # Stop codon
    'TGG': 'W',  # Tryptophan
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',  # Arginine
    'AGT': 'S', 'AGC': 'S',  # Serine
    'AGA': 'R', 'AGG': 'R',  # Arginine
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',  # Glycine
    'REF' : 'REF',
}

# %% Functions
def cmap_red_purple_blues(
    p_red_purple=0.10,      # red to purple limit
    p_transition=0.06,      # purple to blues transition lenght
    left_color="red",
    N=256,
    blues_min=0.2,
    blues_max=1.0
):
    # N colors in each section of the color map
    n_rp = max(2, int(N * p_red_purple))
    n_tr = max(2, int(N * p_transition))
    n_bl = N - n_rp - n_tr

    # Red to purple
    rp_vals = np.linspace(0, 1, n_rp)
    cmap_rp = mcolors.LinearSegmentedColormap.from_list(
        "red_purple", [left_color, "purple"], N=n_rp
    )
    rp_rgba = cmap_rp(rp_vals)

    # Purple to blue soft transition
    blues = cm.get_cmap("Blues")
    first_blue = blues(blues_min)

    tr_vals = np.linspace(0, 1, n_tr)
    cmap_tr = mcolors.LinearSegmentedColormap.from_list(
        "purple_to_blue",
        ["purple", first_blue],
        N=n_tr
    )
    tr_rgba = cmap_tr(tr_vals)

    # Blues original color palette
    bl_rgba = blues(np.linspace(blues_min, blues_max, n_bl))

    # Merged colors to get the final cmap
    colors = np.vstack([
        rp_rgba,
        tr_rgba[1:],   # Avoid doubling first purple
        bl_rgba[1:]    # Avoid doubling first blue
    ])

    return ListedColormap(colors, name="red_purple_Blues")

def aa_label_yaxis(aa_list):
    """
    Given a list like ['A','A','A','R','R','L',...], return
    a list of tuples (aa, start_row, end_row) for contiguous runs.
    """
    runs = []
    if not aa_list:
        return runs
    start = 0
    for i in range(1, len(aa_list) + 1):
        if i == len(aa_list) or aa_list[i] != aa_list[start]:
            runs.append((aa_list[start], start, i - 1))
            start = i
    return runs

# %% Import files
df_ref = pd.read_csv('ref.csv', index_col=0)
gibson_div = pd.read_csv('gibson_barcode_diversity.csv', index_col=0)
gg_div = pd.read_csv('golden_gate_barcode_diversity.csv', index_col=0)
cumul_div = pd.read_csv('cumulative_barcodes_diversity_analysis.csv', index_col=0)

# %% Figure 1B - Cumulative barcode diversity analysis
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (13, 2.66), dpi = 300)
ax1.grid(True, color="lightgray", linewidth=1, alpha=0.8)
ax2.grid(True, color="lightgray", linewidth=1, alpha=0.8)

sns.lineplot(data = cumul_div[cumul_div['Fragment'] == 'F13'], x = 'nclones_bp', marker='o', y = 'perc', ax=ax1, label="% informative barcodes", errorbar="ci")
sns.lineplot(data = cumul_div[cumul_div['Fragment'] == 'F13'], x = 'nclones_bp', marker='o', y = 'perc_mutations', ax=ax1, label = "% mutation coverage", errorbar="ci")
sns.lineplot(data = cumul_div[cumul_div['Fragment'] == 'F13'], x = 'nclones_bp', marker='o', y = '>4', ax=ax1, label = "% barcode diversity for mutations\nrepresented by #barcodes >4", errorbar="ci")
sns.lineplot(data = cumul_div[cumul_div['Fragment'] == 'F13'], x = 'nclones_bp', marker='o', y = '>9', ax=ax1, label = "% barcode diversity for mutations\nrepresented by #barcodes >9", errorbar="ci")

sns.lineplot(data = cumul_div[cumul_div['Fragment'] == 'F43'], x = 'nclones_bp', marker='o', y = 'perc', ax=ax2, label="% informative barcodes", errorbar="ci")
sns.lineplot(data = cumul_div[cumul_div['Fragment'] == 'F43'], x = 'nclones_bp', marker='o', y = 'perc_mutations', ax=ax2, label = "% mutation coverage", errorbar="ci")
sns.lineplot(data = cumul_div[cumul_div['Fragment'] == 'F43'], x = 'nclones_bp', marker='o', y = '>4', ax=ax2, label = "% barcode diversity for mutations\nrepresented by #barcodes >4", errorbar="ci")
sns.lineplot(data = cumul_div[cumul_div['Fragment'] == 'F43'], x = 'nclones_bp', marker='o', y = '>9', ax=ax2, label = "% barcode diversity for mutations\nrepresented by #barcodes >9", errorbar="ci")

ax1.set_title("F13",fontsize=12)
ax1.legend().remove()
ax2.set_title("F43",fontsize=12)

ax1.set_xlabel('Number of transformants per base pair (bp)', fontsize=12)
ax1.set_ylabel('Percentage (%)',fontsize=12)

ax2.set_xlabel('Number of transformants per base pair (bp)',fontsize=12)
ax2.set_ylabel('Percentage (%)',fontsize=12)

ax1.set_ylim(10,100)
ax2.set_ylim(10,100)
ax1.set_xlim(50,1000)
ax2.set_xlim(50,1000)

sns.move_legend(ax2, "upper left", bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig("figure_1B.png", dpi=300, bbox_inches="tight")

# %% Figure 1C - Gibson - Heatmaps nb barcodes per mutated AA unclipped
div_aa = gibson_div[['Fragment', 'ref_aa', 'position', 'mutation_aa','barcode_per_mut_aa']].drop_duplicates()
div_aa_frag = div_aa[(div_aa['Fragment'] == 'F13') | (div_aa['Fragment'] == 'F43')]
fragments = div_aa_frag['Fragment'].unique()
num_fragments = len(fragments)
nplot=1

for start in range(0, num_fragments, nplot):
    subset = fragments[start:start+nplot]
    
    # Define grid size
    cols = math.ceil(math.sqrt(len(subset)))
    rows = math.ceil(len(subset) / cols)

    # Create subplots
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3.8), constrained_layout=True, dpi=1000)
    axes = axes.flatten() if len(subset) > 1 else [axes]

    for i, fragment in enumerate(subset):
        ax = axes[i]
        ax.set_aspect(1)
        
        # Subset data for the fragment
        df_subset = (div_aa_frag[div_aa_frag['Fragment'] == fragment]
                     .pivot(index='mutation_aa', columns='position', values='barcode_per_mut_aa'
        ))
        df_subset = df_subset.loc[aa_order]
        
        # Make color scale
        vmin = 0
        vmax = int(np.nanmax(df_subset.values))
        n_colors = vmax - vmin +1
        red_purple = float((4 + 1)/(vmax))
        transition = float((10 +1 -4)/(vmax))
        cmap = cmap_red_purple_blues(p_red_purple=red_purple, p_transition = transition, N=int(n_colors))
        
        norm = mcolors.BoundaryNorm(
            boundaries=np.arange(vmin, vmax + 2),  # +2 gives correct bin edges
            ncolors=n_colors
        )
        
        # Plot heatmap with fixed color scale
        sns.heatmap(
            df_subset, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
            cbar=True, norm = norm
        )
       
        # Set custom color bar ticks
        cbar = ax.collections[0].colorbar
        ticks = [10] + list(range(0, vmax + 1, 20))
        ticks = [t for t in ticks if vmin <= t <= vmax]
        positions = [t + 0.5 for t in ticks]
        cbar.set_ticks(positions)
        cbar.set_ticklabels(ticks)
        cbar.ax.minorticks_off()
    
        # Titles and labels
        ax.set_title(f"{fragment}")
        ax.set_xlabel("Codon position")
        ax.set_ylabel("")
    
        # ✅ Adjust axis labels
        ax.set_yticks(np.arange(len(df_subset.index)) + 0.5)
        ax.set_yticklabels(df_subset.index, rotation=0, fontproperties=FontProperties(family="monospace"))
    
        xtick_positions = list(range(0, len(df_subset.columns), 14))
        last_pos = int(len(df_subset.columns) - 1)
        if last_pos not in xtick_positions:
            xtick_positions.append(last_pos)  # Append correctly
        xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
        xtick_labels = df_subset.columns[xtick_positions].tolist()  # Ensure a list format
        ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
        ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
        # ✅ Add gray squares for WT aa
        for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
            ref_aa = str(row['ref_aa'])  # Ensure it's a string
            mut_aa = row['position']  # This is a position (number)
        
            # Check if the WT aa (ref_aa) exists in the Y-axis (mutations)
            if ref_aa in df_subset.index and mut_aa in df_subset.columns:
                y_pos = df_subset.index.get_loc(ref_aa)  # Get row index
                x_pos = df_subset.columns.get_loc(mut_aa)  # Get column index
                
                # ✅ Add gray square for WT codons
                ax.add_patch(plt.Rectangle((x_pos + 0.025, y_pos + 0.025), 0.95, 0.95, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Save figures
    fig.savefig(f"figure_1C_{fragment}",bbox_inches="tight", dpi=1000)
    plt.close(fig)

    
# %% Figure 1D - Golden Gate - Heatmaps nb barcodes per mutated AA unclipped
div_aa = gg_div[['Fragment', 'ref_aa', 'position', 'mutation_aa','barcode_per_mut_aa']].drop_duplicates()
div_aa_frag = div_aa[(div_aa['Fragment'] == 'F13') | (div_aa['Fragment'] == 'F43')]
fragments = div_aa_frag['Fragment'].unique()
num_fragments = len(fragments)
nplot=1

for start in range(0, num_fragments, nplot):
    subset = fragments[start:start+nplot]
    
    # Define grid size
    cols = math.ceil(math.sqrt(len(subset)))
    rows = math.ceil(len(subset) / cols)

    # Create subplots
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3.8), constrained_layout=True, dpi=1000)
    axes = axes.flatten() if len(subset) > 1 else [axes]

    for i, fragment in enumerate(subset):
        ax = axes[i]
        ax.set_aspect(1)
        
        # Subset data for the fragment
        df_subset = (div_aa_frag[div_aa_frag['Fragment'] == fragment]
                     .pivot(index='mutation_aa', columns='position', values='barcode_per_mut_aa'
        ))
        df_subset = df_subset.loc[aa_order]
        
        # Make color scale
        vmin = 0
        vmax = int(np.nanmax(df_subset.values))
        n_colors = vmax - vmin +1
        red_purple = float((4 + 1)/(vmax))
        transition = float((10 +1 -4)/(vmax))
        cmap = cmap_red_purple_blues(p_red_purple=red_purple, p_transition = transition, N=int(n_colors))
        
        norm = mcolors.BoundaryNorm(
            boundaries=np.arange(vmin, vmax + 2),  # +2 gives correct bin edges
            ncolors=n_colors
        )
        
        # Plot heatmap with fixed color scale
        sns.heatmap(
            df_subset, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
            cbar=True, norm = norm
        )
       
        # Set custom color bar ticks
        cbar = ax.collections[0].colorbar
        ticks = [10] + list(range(0, vmax + 1, 20))
        ticks = [t for t in ticks if vmin <= t <= vmax]
        positions = [t + 0.5 for t in ticks]
        cbar.set_ticks(positions)
        cbar.set_ticklabels(ticks)
        cbar.ax.minorticks_off()
    
        # Titles and labels
        ax.set_title(f"{fragment}")
        ax.set_xlabel("Codon position")
        ax.set_ylabel("")
    
        # ✅ Adjust axis labels
        ax.set_yticks(np.arange(len(df_subset.index)) + 0.5)
        ax.set_yticklabels(df_subset.index, rotation=0, fontproperties=FontProperties(family="monospace"))
    
        xtick_positions = list(range(0, len(df_subset.columns), 14))
        last_pos = int(len(df_subset.columns) - 1)
        if last_pos not in xtick_positions:
            xtick_positions.append(last_pos)  # Append correctly
        xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
        xtick_labels = df_subset.columns[xtick_positions].tolist()  # Ensure a list format
        ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
        ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
        # ✅ Add gray squares for WT aa
        for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
            ref_aa = str(row['ref_aa'])  # Ensure it's a string
            mut_aa = row['position']  # This is a position (number)
        
            # Check if the WT aa (ref_aa) exists in the Y-axis (mutations)
            if ref_aa in df_subset.index and mut_aa in df_subset.columns:
                y_pos = df_subset.index.get_loc(ref_aa)  # Get row index
                x_pos = df_subset.columns.get_loc(mut_aa)  # Get column index
                
                # ✅ Add gray square for WT codons
                ax.add_patch(plt.Rectangle((x_pos + 0.025, y_pos + 0.025), 0.95, 0.95, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Save figures
    fig.savefig(f"figure_1D_{fragment}",bbox_inches="tight", dpi=1000)
    plt.close(fig)

