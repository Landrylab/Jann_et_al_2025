#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:33:34 2026

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

wkdir = '/home/alicia-pageau/Documents/antifungal_project/PDR1/00_scripts/Jann_et_al_2025/Upadated_scripts_after_review_11_2025/Supplementary_figures/'
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
uniformity_df = pd.read_csv('uniformity.csv', index_col=0)

# %% Supp figure 1 - Gibson - Heatmaps nb barcodes per mutated AA clipped
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
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3.9), constrained_layout=True, dpi=1000)
    axes = axes.flatten() if len(subset) > 1 else [axes]
    
    # Set color scale limits clipped
    vmin, vmax = 1, 10
    n_bins = 10
    tick_step = 1
    
    cmap = mcolors.ListedColormap(plt.cm.Blues(np.linspace(0, 1, n_bins)))
    bounds = np.linspace(vmin, vmax, n_bins + 1)
    tick_positions = (bounds[:-1] + bounds[1:]) / 2
    tick_labels = np.arange(vmin, vmax + tick_step, tick_step)[:len(tick_positions)]
    
    
    for i, fragment in enumerate(subset):
        ax = axes[i]
        ax.set_aspect(1)
        
        # Subset data for the fragment
        df_subset = (div_aa_frag[div_aa_frag['Fragment'] == fragment]
                     .pivot(index='mutation_aa', columns='position', values='barcode_per_mut_aa'
        ))
        df_subset = df_subset.loc[aa_order]
        
        # Get % of clipped barcodes
        barcode_counts = df_subset.fillna(0).values
        total_barcodes = np.sum(barcode_counts)
        clipped_counts = np.sum(np.maximum(barcode_counts - vmax, 0))
        percent_clipped = (clipped_counts / total_barcodes) * 100
        print(f"Percentage of clipped barcodes for {fragment}: {percent_clipped:.2f}%")
        
        # Plot heatmap with fixed color scale
        sns.heatmap(
            df_subset, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
            cbar=True, vmin = vmin, vmax = vmax,
            cbar_kws={"ticks": tick_positions}
        )
        # Set custom tick labels
        cbar = ax.collections[0].colorbar
        cbar.set_ticklabels(tick_labels)
    
        # Titles and labels
        ax.set_title(f"{fragment}")
        ax.set_xlabel("Codon position")
        ax.set_ylabel("")
    
        # Adjust axis labels
        ax.set_yticks(np.arange(len(df_subset.index)) + 0.5)
        ax.set_yticklabels(df_subset.index, rotation=0)
    
        xtick_positions = list(range(0, len(df_subset.columns), 14))
        last_pos = int(len(df_subset.columns) - 1)
        if last_pos not in xtick_positions:
            xtick_positions.append(last_pos)  # Append correctly
        xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
        xtick_labels = df_subset.columns[xtick_positions].tolist()  # Ensure a list format
        ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
        ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
        # Add gray squares for WT aa
        for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
            ref_aa = str(row['ref_aa'])  # Ensure it's a string
            mut_aa = row['position']  # This is a position (number)
        
            # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
            if ref_aa in df_subset.index and mut_aa in df_subset.columns:
                y_pos = df_subset.index.get_loc(ref_aa)  # Get row index
                x_pos = df_subset.columns.get_loc(mut_aa)  # Get column index
                
                # Add gray square for WT codons
                ax.add_patch(plt.Rectangle((x_pos+0.025, y_pos+0.025), 0.95, 0.95, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Save figures
    fig.savefig(f"supp_figure_1_{fragment}",bbox_inches="tight", dpi=1000)
    plt.close(fig)

# %% Supp figure 2 - Golden Gate - Heatmaps nb barcodes per mutated AA clipped
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
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3.9), constrained_layout=True, dpi=1000)
    axes = axes.flatten() if len(subset) > 1 else [axes]
    
    # Set color scale limits clipped
    vmin, vmax = 1, 10
    n_bins = 10
    tick_step = 1
    
    cmap = mcolors.ListedColormap(plt.cm.Blues(np.linspace(0, 1, n_bins)))
    bounds = np.linspace(vmin, vmax, n_bins + 1)
    tick_positions = (bounds[:-1] + bounds[1:]) / 2
    tick_labels = np.arange(vmin, vmax + tick_step, tick_step)[:len(tick_positions)]
    
    
    for i, fragment in enumerate(subset):
        ax = axes[i]
        ax.set_aspect(1)
        
        # Subset data for the fragment
        df_subset = (div_aa_frag[div_aa_frag['Fragment'] == fragment]
                     .pivot(index='mutation_aa', columns='position', values='barcode_per_mut_aa'
        ))
        df_subset = df_subset.loc[aa_order]
        
        # Get % of clipped barcodes
        barcode_counts = df_subset.fillna(0).values
        total_barcodes = np.sum(barcode_counts)
        clipped_counts = np.sum(np.maximum(barcode_counts - vmax, 0))
        percent_clipped = (clipped_counts / total_barcodes) * 100
        print(f"Percentage of clipped barcodes for {fragment}: {percent_clipped:.2f}%")
        
        # Plot heatmap with fixed color scale
        sns.heatmap(
            df_subset, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
            cbar=True, vmin = vmin, vmax = vmax,
            cbar_kws={"ticks": tick_positions}
        )
        # Set custom tick labels
        cbar = ax.collections[0].colorbar
        cbar.set_ticklabels(tick_labels)
    
        # Titles and labels
        ax.set_title(f"{fragment}")
        ax.set_xlabel("Codon position")
        ax.set_ylabel("")
    
        # Adjust axis labels
        ax.set_yticks(np.arange(len(df_subset.index)) + 0.5)
        ax.set_yticklabels(df_subset.index, rotation=0)
    
        xtick_positions = list(range(0, len(df_subset.columns), 14))
        last_pos = int(len(df_subset.columns) - 1)
        if last_pos not in xtick_positions:
            xtick_positions.append(last_pos)  # Append correctly
        xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
        xtick_labels = df_subset.columns[xtick_positions].tolist()  # Ensure a list format
        ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
        ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
        # Add gray squares for WT aa
        for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
            ref_aa = str(row['ref_aa'])  # Ensure it's a string
            mut_aa = row['position']  # This is a position (number)
        
            # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
            if ref_aa in df_subset.index and mut_aa in df_subset.columns:
                y_pos = df_subset.index.get_loc(ref_aa)  # Get row index
                x_pos = df_subset.columns.get_loc(mut_aa)  # Get column index
                
                # Add gray square for WT codons
                ax.add_patch(plt.Rectangle((x_pos+0.025, y_pos+0.025), 0.95, 0.95, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Save figures
    fig.savefig(f"supp_figure_2_{fragment}",bbox_inches="tight", dpi=1000)
    plt.close(fig)
    
# %% Supp figure 3 - Gibson - Heatmaps nb barcodes per mutated codon unclipped
div_frag = gibson_div[(gibson_div['Fragment'] == 'F13') | (gibson_div['Fragment'] == 'F43')]
fragments = div_frag['Fragment'].unique()
num_fragments = len(fragments)
nplot=1

for start in range(0, num_fragments, nplot):
    subset = fragments[start:start+nplot]
    
    # Define grid size
    cols = math.ceil(math.sqrt(len(subset)))
    rows = math.ceil(len(subset) / cols)

    # Create subplots
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5), constrained_layout=True, dpi=1000)
    axes = axes.flatten() if len(subset) > 1 else [axes]
    
    for i, fragment in enumerate(subset):
        ax = axes[i]
        ax.set_aspect(1)
    
        # Subset data for the fragment
        df_subset = div_frag[div_frag['Fragment'] == fragment].pivot(
            index=['mutation'], columns='position', values='barcode_per_mut'
        )
        
        # Order df_subset rows by AA 
        aa_for_row = pd.Series(df_subset.index).map(lambda c: codon_to_aa.get(str(c).upper()))
        order_df = pd.DataFrame({'mutation': df_subset.index, 'aa': aa_for_row.values})
        order_df['aa'] = pd.Categorical(order_df['aa'], categories=aa_order, ordered=True)
        order_df = order_df.sort_values(['aa', 'mutation'], kind='stable')
        aa_labels = order_df['aa'].tolist()
        df_subset = df_subset.loc[order_df['mutation']] 
        
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
    
        # Plot heatmap
        sns.heatmap(
            df_subset, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
            cbar=True, norm = norm
        )
        
        # Set custom color bar ticks
        cbar = ax.collections[0].colorbar
        ticks = [4] + list(range(0, vmax + 1, 10))
        ticks = [t for t in ticks if vmin <= t <= vmax]
        positions = [t + 0.5 for t in ticks]
        cbar.set_ticks(positions)
        cbar.set_ticklabels(ticks)
        cbar.ax.minorticks_off()
    
        # Titles and axis labels
        ax.set_title(f"{fragment}")
        ax.set_xlabel("Mutated Codon")
        ax.set_ylabel("")
    
        # Set custom X axis ticks
        xtick_positions = list(range(0, len(df_subset.columns), 14))
        last_pos = int(len(df_subset.columns) - 1)
        if last_pos not in xtick_positions:
            xtick_positions.append(last_pos)  # Append correctly
        xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
        xtick_labels = df_subset.columns[xtick_positions].tolist()  # Ensure a list format
        ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
        ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
        # Set custom Y axis ticks (to show EVERY CODON)
        n_rows = df_subset.shape[0]
        y_centers = np.arange(n_rows) + 0.5
        ax.yaxis.set_major_locator(FixedLocator(y_centers))  # no decimation
        ax.set_yticklabels(df_subset.index.tolist(), rotation=0,fontsize=8, fontproperties=FontProperties(family="monospace"))
        
        # Add AA label on the left side of the Y axis
        runs = aa_label_yaxis(aa_labels) # Compute AA labels positions to get only one in the center
        x_pad = -3.5  # move left of the y-axis
        for aa, i0, i1 in runs:
            y_mid = (i0 + i1) / 2.0 + 0.5  # center in the block
            ax.text(
                x_pad, y_mid, aa,
                va="center", ha="center",
                fontsize=10, fontweight="bold",
                color="black", clip_on=False,
                fontproperties=FontProperties(family="monospace")
            )
        
        # Add small line to defined each AA block
        x_br = -2.85
        for aa, i0, i1 in runs:
            block_len = (i1 - i0 + 1)
            if block_len <= 1:
                continue  # skip single-codon AA blocks
        
            y0, y1 = i0 + 0.2, i1 + 0.8 
            ax.plot([x_br, x_br], [y0, y1], color="black", lw=1, clip_on=False)
    
        # Add gray squares for WT codons
        for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
            ref_codon = str(row['ref_codon'])  # Ensure it's a string
            mut_codon = row['position']  # This is a position (number)
        
            # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
            if ref_codon in df_subset.index and mut_codon in df_subset.columns:
                y_pos = df_subset.index.get_loc(ref_codon)  # Get row index
                x_pos = df_subset.columns.get_loc(mut_codon)  # Get column index
                
                # ✅ Add gray square for WT codons
                ax.add_patch(plt.Rectangle((x_pos+ 0.025, y_pos + 0.025), 0.95, 0.95, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Save figures
    fig.savefig(f"supp_figure_3_{fragment}",bbox_inches="tight", dpi=1000)
    plt.close(fig)

# %% Supp figure 4 - Golden Gate - Heatmaps nb barcodes per mutated codon unclipped
div_frag = gg_div[(gg_div['Fragment'] == 'F13') | (gg_div['Fragment'] == 'F43')]
fragments = div_frag['Fragment'].unique()
num_fragments = len(fragments)
nplot=1

for start in range(0, num_fragments, nplot):
    subset = fragments[start:start+nplot]
    
    # Define grid size
    cols = math.ceil(math.sqrt(len(subset)))
    rows = math.ceil(len(subset) / cols)

    # Create subplots
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5), constrained_layout=True, dpi=1000)
    axes = axes.flatten() if len(subset) > 1 else [axes]
    
    for i, fragment in enumerate(subset):
        ax = axes[i]
        ax.set_aspect(1)
    
        # Subset data for the fragment
        df_subset = div_frag[div_frag['Fragment'] == fragment].pivot(
            index=['mutation'], columns='position', values='barcode_per_mut'
        )
        
        # Order df_subset rows by AA 
        aa_for_row = pd.Series(df_subset.index).map(lambda c: codon_to_aa.get(str(c).upper()))
        order_df = pd.DataFrame({'mutation': df_subset.index, 'aa': aa_for_row.values})
        order_df['aa'] = pd.Categorical(order_df['aa'], categories=aa_order, ordered=True)
        order_df = order_df.sort_values(['aa', 'mutation'], kind='stable')
        aa_labels = order_df['aa'].tolist()
        df_subset = df_subset.loc[order_df['mutation']] 
        
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
    
        # Plot heatmap
        sns.heatmap(
            df_subset, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
            cbar=True, norm = norm
        )
        
        # Set custom color bar ticks
        cbar = ax.collections[0].colorbar
        ticks = [4] + list(range(0, vmax + 1, 10))
        ticks = [t for t in ticks if vmin <= t <= vmax]
        positions = [t + 0.5 for t in ticks]
        cbar.set_ticks(positions)
        cbar.set_ticklabels(ticks)
        cbar.ax.minorticks_off()
    
        # Titles and axis labels
        ax.set_title(f"{fragment}")
        ax.set_xlabel("Mutated Codon")
        ax.set_ylabel("")
    
        # Set custom X axis ticks
        xtick_positions = list(range(0, len(df_subset.columns), 14))
        last_pos = int(len(df_subset.columns) - 1)
        if last_pos not in xtick_positions:
            xtick_positions.append(last_pos)  # Append correctly
        xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
        xtick_labels = df_subset.columns[xtick_positions].tolist()  # Ensure a list format
        ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
        ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
        # Set custom Y axis ticks (to show EVERY CODON)
        n_rows = df_subset.shape[0]
        y_centers = np.arange(n_rows) + 0.5
        ax.yaxis.set_major_locator(FixedLocator(y_centers))  # no decimation
        ax.set_yticklabels(df_subset.index.tolist(), rotation=0,fontsize=8, fontproperties=FontProperties(family="monospace"))
        
        # Add AA label on the left side of the Y axis
        runs = aa_label_yaxis(aa_labels) # Compute AA labels positions to get only one in the center
        x_pad = -3.5  # move left of the y-axis
        for aa, i0, i1 in runs:
            y_mid = (i0 + i1) / 2.0 + 0.5  # center in the block
            ax.text(
                x_pad, y_mid, aa,
                va="center", ha="center",
                fontsize=10, fontweight="bold",
                color="black", clip_on=False,
                fontproperties=FontProperties(family="monospace")
            )
        
        # Add small line to defined each AA block
        x_br = -2.85
        for aa, i0, i1 in runs:
            block_len = (i1 - i0 + 1)
            if block_len <= 1:
                continue  # skip single-codon AA blocks
        
            y0, y1 = i0 + 0.2, i1 + 0.8 
            ax.plot([x_br, x_br], [y0, y1], color="black", lw=1, clip_on=False)
    
        # Add gray squares for WT codons
        for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
            ref_codon = str(row['ref_codon'])  # Ensure it's a string
            mut_codon = row['position']  # This is a position (number)
        
            # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
            if ref_codon in df_subset.index and mut_codon in df_subset.columns:
                y_pos = df_subset.index.get_loc(ref_codon)  # Get row index
                x_pos = df_subset.columns.get_loc(mut_codon)  # Get column index
                
                # ✅ Add gray square for WT codons
                ax.add_patch(plt.Rectangle((x_pos+ 0.025, y_pos + 0.025), 0.95, 0.95, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Save figures
    fig.savefig(f"supp_figure_4_{fragment}",bbox_inches="tight", dpi=1000)
    plt.close(fig)

# %% Supp figure 5 - Gibson - Heatmap nb barcodes per mutated AA unclipped whole PDR1 sequence
div_aa = gibson_div[['Fragment', 'ref_aa', 'position', 'mutation_aa','barcode_per_mut_aa']].drop_duplicates()
df = div_aa.pivot(index='mutation_aa', columns='position', values='barcode_per_mut_aa')
df = df.loc[aa_order]

fig, ax = plt.subplots(figsize=(200, 10), dpi=225)
# Make color scale
vmin = 0
vmax = int(np.nanmax(df.values))
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
    df, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
    cbar=True, vmin = vmin, vmax = vmax,
    cbar_kws={"pad": 0.005}
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
ax.set_xlabel("Codon position")
ax.set_ylabel("")

# ✅ Adjust axis labels
ax.set_yticks(np.arange(len(df.index)) + 0.5)
ax.set_yticklabels(df.index, rotation=0, fontproperties=FontProperties(family="monospace"))

xtick_positions = list(range(0, len(df.columns), 25))
xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
xtick_labels = df.columns[xtick_positions].tolist()  # Ensure a list format
ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
# ✅ Add gray squares for WT codons
for idx, row in df_ref.iterrows():
    ref_aa = str(row['ref_aa'])  # Ensure it's a string
    mut_aa = row['position']  # This is a position (number)
        
    # Check if the WT aa (ref_aa) exists in the Y-axis (mutations)
    if ref_aa in df.index and mut_aa in df.columns:
        y_pos = df.index.get_loc(ref_aa)  # Get row index
        x_pos = df.columns.get_loc(mut_aa)  # Get column index
                
        # ✅ Add gray square for WT aa
        ax.add_patch(plt.Rectangle((x_pos + 0.025, y_pos + 0.015), 0.95, 0.97, color='gray', lw=0.1))
    
# Save plot
plt.savefig("supp_figure_5.png", dpi=225, bbox_inches="tight")

# %% Supp figure 6 - Golden Gate - Heatmap nb barcodes per mutated AA unclipped whole PDR1 sequence
div_aa = gg_div[['Fragment', 'ref_aa', 'position', 'mutation_aa','barcode_per_mut_aa']].drop_duplicates()
df = div_aa.pivot(index='mutation_aa', columns='position', values='barcode_per_mut_aa')
df = df.loc[aa_order]

fig, ax = plt.subplots(figsize=(200, 10), dpi=225)
# Make color scale
vmin = 0
vmax = int(np.nanmax(df.values))
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
    df, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
    cbar=True, vmin = vmin, vmax = vmax,
    cbar_kws={"pad": 0.005}
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
ax.set_xlabel("Codon position")
ax.set_ylabel("")

# ✅ Adjust axis labels
ax.set_yticks(np.arange(len(df.index)) + 0.5)
ax.set_yticklabels(df.index, rotation=0, fontproperties=FontProperties(family="monospace"))

xtick_positions = list(range(0, len(df.columns), 25))
xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
xtick_labels = df.columns[xtick_positions].tolist()  # Ensure a list format
ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
# ✅ Add gray squares for WT codons
for idx, row in df_ref.iterrows():
    ref_aa = str(row['ref_aa'])  # Ensure it's a string
    mut_aa = row['position']  # This is a position (number)
        
    # Check if the WT aa (ref_aa) exists in the Y-axis (mutations)
    if ref_aa in df.index and mut_aa in df.columns:
        y_pos = df.index.get_loc(ref_aa)  # Get row index
        x_pos = df.columns.get_loc(mut_aa)  # Get column index
                
        # ✅ Add gray square for WT aa
        ax.add_patch(plt.Rectangle((x_pos + 0.025, y_pos + 0.015), 0.95, 0.97, color='gray', lw=0.1))
    
# Save plot
plt.savefig("supp_figure_6.png", dpi=225, bbox_inches="tight")

# %% Supp figure 7 - Gibson - Heatmap nb barcodes per mutated codon unclipped whole PDR1 sequence
df = gibson_div.pivot(index='mutation', columns='position', values='barcode_per_mut')
df = df[~df.index.isna()]
df = df.fillna(0)

fig, ax = plt.subplots(figsize=(200, 10), dpi=225)

# Order df rows by AA 
aa_for_row = pd.Series(df.index).map(lambda c: codon_to_aa.get(str(c).upper()))
order_df = pd.DataFrame({'mutation': df.index, 'aa': aa_for_row.values})
order_df['aa'] = pd.Categorical(order_df['aa'], categories=aa_order, ordered=True)
order_df = order_df.sort_values(['aa', 'mutation'], kind='stable')
aa_labels = order_df['aa'].tolist()
df = df.loc[order_df['mutation']] 

# Make color scale
vmin = 0
vmax = int(np.nanmax(df.values))
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
    df, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
    cbar=True, vmin = vmin, vmax = vmax,
    cbar_kws={"pad": 0.005}
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
ax.set_xlabel("Codon position")
ax.set_ylabel("")

# ✅ Adjust axis labels
ax.set_yticks(np.arange(len(df.index)) + 0.5)
ax.set_yticklabels(df.index, rotation=0, fontproperties=FontProperties(family="monospace"))

xtick_positions = list(range(0, len(df.columns), 25))
xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
xtick_labels = df.columns[xtick_positions].tolist()  # Ensure a list format
ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
# ✅ Add gray squares for WT codons
for idx, row in df_ref.iterrows():
    ref_codon = str(row['ref_codon'])  # Ensure it's a string
    mut_codon = row['position']  # This is a position (number)
        
    # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
    if ref_codon in df.index and mut_codon in df.columns:
        y_pos = df.index.get_loc(ref_codon)  # Get row index
        x_pos = df.columns.get_loc(mut_codon)  # Get column index
                
        # ✅ Add gray square for WT codons
        ax.add_patch(plt.Rectangle((x_pos + 0.025, y_pos + 0.015), 0.95, 0.97, color='gray', lw=0.1))

# Save plot
plt.savefig("supp_figure_7.png", dpi=225, bbox_inches="tight")

# %% Supp figure 8 - Golden Gate - Heatmap nb barcodes per mutated codon unclipped whole PDR1 sequence
df = gg_div.pivot(index='mutation', columns='position', values='barcode_per_mut')
df = df[~df.index.isna()]
df = df.fillna(0)

fig, ax = plt.subplots(figsize=(200, 10), dpi=225)

# Order df rows by AA 
aa_for_row = pd.Series(df.index).map(lambda c: codon_to_aa.get(str(c).upper()))
order_df = pd.DataFrame({'mutation': df.index, 'aa': aa_for_row.values})
order_df['aa'] = pd.Categorical(order_df['aa'], categories=aa_order, ordered=True)
order_df = order_df.sort_values(['aa', 'mutation'], kind='stable')
aa_labels = order_df['aa'].tolist()
df = df.loc[order_df['mutation']] 

# Make color scale
vmin = 0
vmax = int(np.nanmax(df.values))
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
    df, cmap=cmap, annot=False, linewidths=0.5, ax=ax,
    cbar=True, vmin = vmin, vmax = vmax,
    cbar_kws={"pad": 0.005}
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
ax.set_xlabel("Codon position")
ax.set_ylabel("")

# ✅ Adjust axis labels
ax.set_yticks(np.arange(len(df.index)) + 0.5)
ax.set_yticklabels(df.index, rotation=0, fontproperties=FontProperties(family="monospace"))

xtick_positions = list(range(0, len(df.columns), 25))
xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
xtick_labels = df.columns[xtick_positions].tolist()  # Ensure a list format
ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
# ✅ Add gray squares for WT codons
for idx, row in df_ref.iterrows():
    ref_codon = str(row['ref_codon'])  # Ensure it's a string
    mut_codon = row['position']  # This is a position (number)
        
    # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
    if ref_codon in df.index and mut_codon in df.columns:
        y_pos = df.index.get_loc(ref_codon)  # Get row index
        x_pos = df.columns.get_loc(mut_codon)  # Get column index
                
        # ✅ Add gray square for WT codons
        ax.add_patch(plt.Rectangle((x_pos + 0.025, y_pos + 0.015), 0.95, 0.97, color='gray', lw=0.1))

# Save plot
plt.savefig("supp_figure_8.png", dpi=225, bbox_inches="tight")

# %% Supp figure 9 - Gini boxplot
plt.figure(figsize=(3, 6), dpi=300)
sns.boxplot(data=uniformity_df, x='exp', y='Gini',fill=False, 
            color='black',fliersize=0, order=['gibson','golden gate'])
sns.stripplot(data=uniformity_df, x='exp', y='Gini', 
              marker='o', size=6, linewidth=0.5, color='black', 
              alpha=0.6, order=['gibson','golden gate'])
plt.ylim(0,1)
plt.xticks(ticks=[0, 1], labels=["Gibson", "Golden Gate"])
plt.ylabel('Uniformity (Gini coefficient)')
plt.xlabel('')
plt.savefig("supp_figure_9.png", dpi=300, bbox_inches="tight")


# %% Supp figure 10 - Uniformity boxplot
plt.figure(figsize=(3, 6), dpi=300)
sns.boxplot(data=uniformity_df, x='exp', y='LogDiff',fill=False, 
            color='black',fliersize=0, order=['gibson','golden gate'])
sns.stripplot(data=uniformity_df, x='exp', y='LogDiff', 
              marker='o', size=6, linewidth=0.5, color='black', 
              alpha=0.6, order=['gibson','golden gate'])
plt.ylim(0)
plt.xticks(ticks=[0, 1], labels=["Gibson", "Golden Gate"])
plt.ylabel('Uniformity (LogDiff)')
plt.xlabel('')
plt.savefig("supp_figure_10.png", dpi=300, bbox_inches="tight")

# %% Supp figure 11 - Distribution nb barcodes per mutated AA
# Gibson
div_aa = gibson_div[['Fragment', 'ref_aa', 'position', 'mutation_aa','barcode_per_mut_aa']].drop_duplicates()
div_aa_frag = div_aa[(div_aa['Fragment'] == 'F13') | (div_aa['Fragment'] == 'F43')]

palette = {
    'F13': '#666666',  # gray
    'F43': '#56B4E9',  # blue
}

plt.figure(figsize=(10, 6))
sns.set(rc={'axes.facecolor': 'white',
 'axes.edgecolor': 'black',
 'axes.grid': True,
 'figure.facecolor': 'white',
 'grid.color': '#b0b0b0',
 'xtick.direction': 'out',
 'ytick.direction': 'out',
 'xtick.bottom': True,
 'ytick.left': True,
 })

sns.histplot(
    data=div_aa_frag,
    x="barcode_per_mut_aa",
    hue="Fragment",
    binwidth=2,
    alpha=0.5,
    palette=palette
)

plt.xticks(range(0, int(div_aa_frag["barcode_per_mut_aa"].max()) + 10, 10))
plt.xlabel("Number of barcode per mutation (AA)")
plt.ylabel("Count")
plt.savefig("supp_figure_11_gibson.png", dpi=300, bbox_inches="tight")

# Golden Gate
div_aa = gg_div[['Fragment', 'ref_aa', 'position', 'mutation_aa','barcode_per_mut_aa']].drop_duplicates()
div_aa_frag = div_aa[(div_aa['Fragment'] == 'F13') | (div_aa['Fragment'] == 'F43')]

palette = {
    'F13': '#666666',  # gray
    'F43': '#56B4E9',  # blue
}

plt.figure(figsize=(10, 6))
sns.set(rc={'axes.facecolor': 'white',
 'axes.edgecolor': 'black',
 'axes.grid': True,
 'figure.facecolor': 'white',
 'grid.color': '#b0b0b0',
 'xtick.direction': 'out',
 'ytick.direction': 'out',
 'xtick.bottom': True,
 'ytick.left': True,
 })

sns.histplot(
    data=div_aa_frag,
    x="barcode_per_mut_aa",
    hue="Fragment",
    binwidth=2,
    alpha=0.5,
    palette=palette
)

plt.xticks(range(0, int(div_aa_frag["barcode_per_mut_aa"].max()) + 10, 10))
plt.xlabel("Number of barcode per mutation (AA)")
plt.ylabel("Count")
plt.savefig("supp_figure_11_golden_gate.png", dpi=300, bbox_inches="tight")