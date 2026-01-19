#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 09:56:26 2025

10-11-2025
PDR1 Golden Gate Sequencing Analysis
"""

import pandas as pd
import numpy as np
import seaborn as sns
import Bio
from Bio import SeqIO
import matplotlib.pyplot as plt
import glob
import os
import re
import math
import natsort
from natsort import natsorted
from matplotlib import colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FixedLocator
from matplotlib.font_manager import FontProperties
import gzip

print('pandas', pd.__version__)
print('matplotlib', plt.matplotlib.__version__)
print('numpy', np.__version__)
print('seaborn', sns.__version__)
print('Bio', Bio.__version__)
print('re', re.__version__)
print('natsort', natsort.__version__)

wkdir = '/home/alicia-pageau/Documents/antifungal_project/PDR1'
os.chdir(wkdir)

#%% Functions

def check_codons(seq1, seq2):
    n_codons = np.floor(len(seq1)/3)
    changed_codons = []
    referen_codons = []
    for codon in range(0, len(seq1), 3):
        codon_seq1 = seq1[codon:codon+3]
        codon_seq2 = seq2[codon:codon+3]
        diffs = [i for i,x in enumerate(zip(codon_seq1,codon_seq2)) if x[0]!=x[1]]
        n_diffs = len(diffs)
        if n_diffs > 0:
            changed_codons.append(codon_seq1)
            referen_codons.append(codon_seq2)
    if changed_codons:
        codon_list = ','.join(changed_codons)
        ref_list = ','.join(referen_codons)
    else:
        codon_list = 'REF'
        ref_list = 'REF'
    codon_number = len(changed_codons)
    return(ref_list, codon_list, codon_number)

def mutation_positions(sequence1, sequence2):
    pos = [i+1 for i,x in enumerate(zip(sequence1,sequence2)) if x[0]!=x[1]]
    if pos:
        return pos
    else:
        return [0]

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

# aa order per properties for heatmaps
aa_order = ["*", "P", "G", "C", "Q", "N", "T", "S", "E", "D",
            "K", "H", "R", "W", "Y", "F", "M", "L", "I", "V", "A"]

def is_valid_nnk(codon):
    return bool(re.match(r'^[ACGT]{2}[GT]$', codon))  # First two bases: A/C/G/T, Last base: G/T

def gini(array):
    """Compute Gini coefficient of array of counts, treating NaN as 0."""
    array = np.nan_to_num(array, nan=0.0)
    if np.all(array == 0):
        return np.nan
    sorted_array = np.sort(array)
    n = len(array)
    cumvals = np.cumsum(sorted_array)
    return (n + 1 - 2 * np.sum(cumvals) / cumvals[-1]) / n

def compute_uniformity(df, count_col='reads', group_col='Fragment'):
    """Compute CV, Shannon entropy, Gini, and LogDiff per group (fragment), treating NaN as 0."""
    results = []
    for name, group in df.groupby(group_col):
        counts = group[count_col].values
        counts = np.nan_to_num(counts, nan=0.0)  # replace NaN with 0

        mean = counts.mean()
        std = counts.std()
        cv = std / mean if mean > 0 else np.nan

        total_counts = counts.sum()
        freq = counts / total_counts if total_counts > 0 else np.zeros_like(counts)

        # Shannon uniformity
        entropy = -(freq * np.log2(freq + 1e-12)).sum()  # avoid log(0)
        max_entropy = np.log2(len(freq)) if len(freq) > 0 else np.nan
        uniformity = entropy / max_entropy if max_entropy > 0 else np.nan

        # Gini
        g = gini(counts)

        # Log difference between 90th and 10th percentiles
        if len(counts) > 1:
            p90 = np.percentile(counts, 90)
            p10 = np.percentile(counts, 10)
            if p90 > 0 and p10 > 0:
                logdiff = np.log10(p90) - np.log10(p10)
            else:
                logdiff = np.nan
        else:
            logdiff = np.nan

        results.append({
            group_col: name,
            'CV': cv,
            'Shannon_uniformity': uniformity,
            'Gini': g,
            'LogDiff': logdiff
        })

    return pd.DataFrame(results)

# Function for plotting
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
#%% Reference (expected) barcode
df = pd.read_excel("pdr1_barcode_combinations_full_corrected.xlsx")

df_len = df[['Fragment','tot_len','extra_front','extra_back','frag_len']].reset_index(drop=True)
barcode_length = 30
tot_len = df_len.set_index('Fragment')['tot_len'].to_dict()
frag_len = df_len.set_index('Fragment')['frag_len'].to_dict()
extra_front = df_len.set_index('Fragment')['extra_front'].to_dict()
extra_back = df_len.set_index('Fragment')['extra_back'].to_dict()

dbarc = df[['Fragment','codon','pattern']].reset_index(drop=True)

dref = list(SeqIO.parse("references_PDR1_corrected.fasta","fasta"))

dpat = df[['Fragment','pattern']]

R = {rec.id.split('_')[1]: str(rec.seq).upper() for rec in dref}
D = {rec.id.split('_')[1]:len(rec.seq) for rec in dref}

aa_data = []
for seq_id, sequence in R.items():
    for i in range(0, len(sequence) - 2, 3):  # Iterate in steps of 3
        codon = sequence[i:i+3]  # Extract codon
        codon_position = i // 3 + 1  # Convert base position to codon number
        
        # 🔹 Convert codon to amino acid 
        amino_acid = codon_to_aa.get(codon)
        
        aa_data.append([seq_id, codon_position, codon, amino_acid])

# Create DataFrame
df_ref = pd.DataFrame(aa_data, columns=["Fragment", "mut_codon", "ref_codon", "ref_aa"])
df_ref['position'] = ((df_ref['Fragment'].str.replace("F","").astype(int) -1)*25) + df_ref['mut_codon']

#%% Processing reads
reads = glob.glob(f"{wkdir}/golden_gate/04_merged/Pool*.fasta.gz") + \
    glob.glob(f"{wkdir}/golden_gate/04_merged/01_2025/GoldenGate*.fasta.gz") + \
    glob.glob(f"{wkdir}/golden_gate/04_merged/09_2025/GG_*.fasta.gz")

STATS = {}
T1 = []
T2 = []
M = []

for fread in reads: 
    print(fread)
    #fragment = "F" + fread.split('/')[-1].split('.')[0].replace("GoldenGate","")
    fragment = fread.split('/')[-1].split('_')[1]

    open_func = gzip.open if fread.endswith(".gz") else open
    with open_func(fread, "rt") as handle:
       dread = [str(rec.seq) for rec in SeqIO.parse(handle, "fasta")]

    # Reading and processing of reads
    dR = pd.DataFrame({"barcode" : dread})
    a = dR.shape[0]
    print(dR.shape)
    if dR.shape[0] > 0:
        # Removing unaligned
        dR['len'] = dR['barcode'].apply(lambda x: len(x))
        dR = dR[dR['len'] == barcode_length].reset_index(drop=True)
        b = dR.shape[0]
        print(dR.shape)
        
        if dR.shape[0] > 0:
            # Subset pattern to match with ref and remove reads with mismatches in barcode pattern
            dR['pattern'] = dR['barcode'].apply(lambda x: x[2:14] + "_" + x[16:28])
            dR['Fragment'] = fragment 
            dR['mismatch'] = dR.merge(dpat, on=['pattern', 'Fragment'], how='left', indicator=True)['_merge'] == 'both'
            dR = dR[dR['mismatch'] == True]
            c = dR.shape[0]
            print(dR.shape)
            if dR.shape[0] > 0:
                # Removing sequences with N in the sequence
                dR = dR[~dR['barcode'].str.contains('N')].reset_index(drop=True)
                d = dR.shape[0]
                print(dR.shape)
                
                # Count reads per barcode (unique barcode could be sequencing error)
                dR = dR.groupby(['Fragment','barcode','pattern']).agg(barcode_reads = ('barcode','count')).reset_index()
                dR['barcode_reads_>1'] = dR['barcode_reads'] > 1
                dR['barcode_reads_>=4'] = dR['barcode_reads'] >= 4
                dR['barcode_reads_>=9'] = dR['barcode_reads'] >= 9
                e = dR['barcode_reads_>1'].sum()
                f = dR['barcode_reads_>=4'].sum()
                g = dR['barcode_reads_>=9'].sum()
                print(dR.shape)

                if dR.shape[0] > 0:
                    # Merging reads with barcode patterns
                    dM = dR.copy().drop(columns={'barcode_reads_>1','barcode_reads_>=4','barcode_reads_>=9'})
                    dM.to_csv(f"{wkdir}/golden_gate/05_analyze/09_2025/"+fragment+".csv", sep=",", header=True, index=False)
                    M.append(dM)
                    
                    print(dM.shape)
                    STATS[fragment] = [a,b,c,d,e,f,g]

dSTAT = pd.DataFrame(STATS).T.reset_index()
dSTAT.columns = ['Fragment', 'N_reads','N_aligned','N_barcode','N_noN','N_noUnique','N_>=4','N_>=9']
dSTAT['retained'] = dSTAT.apply(lambda x: (x['N_noUnique']/x['N_reads'])*100, axis=1)
dSTAT = dSTAT.sort_values(by=['Fragment']).reset_index(drop=True)
dSTAT.to_csv(f"{wkdir}/golden_gate/05_analyze/09_2025/goldenGate_F38_V2_STATS.tsv", sep='\t', header=True, index=False)


#%% Merge with barcode-mutation association from gibson analysis
read_depth_threshold = 1

gibson = pd.read_csv(f"{wkdir}/gibson_all_fragments_sequencing_read_20102025.csv")

golden_gate_F1_F13_F43 = pd.read_csv(f"{wkdir}/golden_gate/05_analyze/GoldenGate_sequencing_reads.csv") 
golden_gate_01_2025 = pd.read_csv(f"{wkdir}/golden_gate/05_analyze/01_2025/GoldenGate_sequencing_reads.csv")
golden_gate_01_2025['GG_assembly'] = 1
golden_gate_01_2025['Pool'] = 1
golden_gate_09_2025 = pd.read_csv(f"{wkdir}/golden_gate/05_analyze/09_2025/GoldenGate_sequencing_reads_with_F38_v2.csv")
golden_gate_09_2025['GG_assembly'] = 2
golden_gate_09_2025['Pool'] = 1
golden_gate_09_2025.loc[golden_gate_09_2025['Fragment'] == 'F1', 'GG_assembly'] = 16
golden_gate_09_2025.loc[golden_gate_09_2025['Fragment'] == 'F1', 'Pool'] = golden_gate_09_2025.loc[golden_gate_09_2025['Fragment'] == 'F1', 'Pool'].fillna(1)

golden_gate = pd.concat([golden_gate_F1_F13_F43,golden_gate_01_2025,golden_gate_09_2025], ignore_index=True)
golden_gate = golden_gate[golden_gate['barcode_reads'] > read_depth_threshold]


# Assign 'pcr' values
pcr_mapping = pd.read_csv(f"{wkdir}/GG_gisbon_pcr_mapping.csv", index_col=0)
golden_gate = golden_gate.merge(pcr_mapping, on=["Fragment", "GG_assembly", "Pool"], how="left")
data = pd.merge(gibson, golden_gate.drop(columns='pattern'), how='inner')

# Barcode uniqueness - Get barcode associated to only one mutation
data['barcode_uniqueness'] = data.groupby(['Fragment','barcode'], dropna=False)['core'].transform('nunique')
data_informative = data[data['barcode_uniqueness'] == 1].drop(columns='barcode_uniqueness').drop_duplicates()
data_informative = data_informative.drop(columns='seq').drop_duplicates()

# Sum reads the different golden gate of the same fragment
cols = ['len', 'core', 'barcode', 'pattern', 'ref', 'mismatches',
       'ref_codon', 'mutation', 'changed_codons_n', 'mut_pos', 'mut_codon',
       'Fragment', 'codon']
data_informative_sum = data_informative.groupby(cols,dropna=False).agg(barcode_reads = ('barcode_reads','sum')).reset_index()
data_informative_sum['barcode_reads_>1'] = data_informative_sum['barcode_reads'] > 1
data_informative_sum['barcode_reads_>=4'] = data_informative_sum['barcode_reads'] >= 4
data_informative_sum['barcode_reads_>=9'] = data_informative_sum['barcode_reads'] >= 9

data_informative_sum['mutation_aa'] = data_informative_sum['mutation'].map(codon_to_aa)
data_informative_sum['ref_aa'] = data_informative_sum['ref_codon'].map(codon_to_aa)

nb_informative_barcode = (
    data_informative_sum
    .groupby('Fragment')
    .agg(nb_informative_barcode=('barcode', 'nunique')).reset_index()
)
nb_informative_barcode.to_csv('/home/alicia-pageau/Documents/antifungal_project/PDR1/golden_gate/05_analyze/golden_gate_informative_barcodes.csv')

#%% Figures (codon)
barcode_mut = data_informative_sum.groupby(['Fragment','mutation','mut_codon']).agg(barcode_per_mut = ('barcode','nunique'), reads = ('barcode_reads','sum')).reset_index()
barcode_mut = barcode_mut[barcode_mut['mutation'] != 'REF']
barcode_mut = barcode_mut[barcode_mut['mutation'].apply(is_valid_nnk)]

# Add position mutation in the whole protein seq
barcode_mut['position'] = ((barcode_mut['Fragment'].str.replace("F","").astype(int) -1)*25) + barcode_mut['mut_codon']

# Create a full DataFrame with all combinations
fragments = set(barcode_mut['Fragment'].unique())  # Convert to set for faster lookup
full_index = [(frag, codon) for frag, length in frag_len.items() if frag in fragments for codon in range(1, (length//3) + 1)]
full_df = pd.DataFrame(full_index, columns=['Fragment', 'mut_codon'])
full_df = pd.merge(full_df, df_ref, how = 'left')

barcode_mut = full_df.merge(barcode_mut, how='left')
barcode_mut['Fragment'] = pd.Categorical(barcode_mut['Fragment'], categories=natsorted(barcode_mut['Fragment'].unique()), ordered=True)
barcode_mut = barcode_mut.sort_values(['Fragment', 'mut_codon'])
barcode_mut['mutation_aa'] = barcode_mut['mutation'].map(codon_to_aa)
barcode_mut['ref_aa'] = barcode_mut['ref_codon'].map(codon_to_aa)

# # In prevision of plot with aa as well
# aa_counts = (
#     barcode_mut
#     .groupby(['position', 'mutation_aa'])
#     .agg(barcode_per_mut_aa=('barcode_per_mut', 'sum'))
#     .reset_index()
# )
# barcode_mut = barcode_mut.merge(aa_counts, on=['position', 'mutation_aa'], how='left')
# barcode_mut.to_csv('/home/alicia-pageau/Documents/antifungal_project/PDR1/00_scripts/Jann_et_al_2025/Upadated_scripts_after_review_11_2025/Figure1/golden_gate_barcode_diversity.csv')

## Compute an plot uniformity
uniformity_df = compute_uniformity(barcode_mut)
#uniformity_df.to_csv('/home/alicia-pageau/Documents/antifungal_project/PDR1/golden_gate/05_analyze/golden_gate_uniformity.csv')
uniformity_df['exp'] = 'golden gate'
gibson_uniformity_df = pd.read_csv(f"{wkdir}/gibson/05_analyze/gibson_uniformity.csv")
gibson_uniformity_df['exp'] = 'gibson'
uniformity_df = pd.concat([uniformity_df,gibson_uniformity_df])

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
plt.show()

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
plt.show()

## PLOT nb barcode per codon whole seq
fig, ax = plt.subplots(figsize=(200, 10), dpi=225)

df = barcode_mut.pivot(index='mutation', columns='position', values='barcode_per_mut')
df = df.fillna(0)
df = df[~df.index.isna()]

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

# Show plot
plt.show()

## PLOT nb barcode per codon
barcode_mut = barcode_mut[(barcode_mut['Fragment'] == 'F13') | (barcode_mut['Fragment'] == 'F43')]
fragments = barcode_mut['Fragment'].unique()
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
        df_subset = barcode_mut[barcode_mut['Fragment'] == fragment].pivot(
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
        ax.set_yticklabels(df_subset.index.tolist(), rotation=0, fontsize=6, fontproperties=FontProperties(family="monospace"))
        
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
                ax.add_patch(plt.Rectangle((x_pos+ 0.05, y_pos + 0.05), 0.90, 0.90, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Show plots
    plt.show()

## PLOT distribution nb barcode per codon
barcode_mut = barcode_mut[(barcode_mut['Fragment'] == 'F13') | (barcode_mut['Fragment'] == 'F43')]
barcode_mut['Fragment'] = barcode_mut['Fragment'].cat.remove_unused_categories()

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
    data=barcode_mut,
    x="barcode_per_mut",
    hue="Fragment",
    binwidth=1,
    alpha=0.5,
    palette=palette
)

plt.xticks(range(0, int(barcode_mut["barcode_per_mut"].max()) + 10, 10))
plt.xlabel("Number of barcode per mutation")
plt.ylabel("Count")
plt.show()

#%% Figures (AA)
barcode_mut = data_informative_sum.groupby(['Fragment','mutation_aa','mut_codon']).agg(barcode_per_mut = ('barcode','nunique'), reads = ('barcode_reads','sum')).reset_index()
barcode_mut = barcode_mut[barcode_mut['mutation_aa'] != 'REF']
barcode_mut = barcode_mut.fillna(0)

plt.figure(figsize=(6, 4), dpi=300)
sns.scatterplot(data=barcode_mut, x="barcode_per_mut", y = "reads", hue='Fragment',legend=None)
plt.xlabel("Number of barcodes per mutation")
plt.ylabel("Number of reads per mutation")
plt.show()

# Add position mutation in the whole protein seq
barcode_mut['position'] = ((barcode_mut['Fragment'].str.replace("F","").astype(int) -1)*25) + barcode_mut['mut_codon']

# Create a full DataFrame with all combinations
fragments = set(barcode_mut['Fragment'].unique())  # Convert to set for faster lookup
full_index = [(frag, codon) for frag, length in frag_len.items() if frag in fragments for codon in range(1, (length//3) + 1)]
full_df = pd.DataFrame(full_index, columns=['Fragment', 'mut_codon'])
full_df = pd.merge(full_df, df_ref, how = 'left')

barcode_mut = full_df.merge(barcode_mut, how='left')
barcode_mut['Fragment'] = pd.Categorical(barcode_mut['Fragment'], categories=natsorted(barcode_mut['Fragment'].unique()), ordered=True)
barcode_mut = barcode_mut.sort_values(['Fragment', 'mut_codon'])

## PLOT nb barcode per mutation whole seq
fig, ax = plt.subplots(figsize=(200, 10), dpi=225)

df = barcode_mut.pivot(index='mutation_aa', columns='position', values='barcode_per_mut')
df = df.loc[aa_order]

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
ax.set_yticklabels(df.index, rotation=0)

xtick_positions = list(range(0, len(df.columns), 25))
xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
xtick_labels = df.columns[xtick_positions].tolist()  # Ensure a list format
ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
        
# ✅ Add gray squares for WT aa
for idx, row in df_ref.iterrows():
    ref_aa = str(row['ref_aa'])  # Ensure it's a string
    mut_aa = row['position']  # This is a position (number)
        
    # Check if the WT codon (ref_aa) exists in the Y-axis (mutations)
    if ref_aa in df.index and mut_aa in df.columns:
        y_pos = df.index.get_loc(ref_aa)  # Get row index
        x_pos = df.columns.get_loc(mut_aa)  # Get column index
                
        # ✅ Add gray square for WT aa
        ax.add_patch(plt.Rectangle((x_pos+ 0.025, y_pos + 0.015), 0.95, 0.97, color='gray', lw=0.1))
    
# Show plots
plt.show()

## PLOT distribution nb barcode per codon
barcode_mut = barcode_mut[(barcode_mut['Fragment'] == 'F13') | (barcode_mut['Fragment'] == 'F43')]
barcode_mut['Fragment'] = barcode_mut['Fragment'].cat.remove_unused_categories()

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
    data=barcode_mut,
    x="barcode_per_mut",
    hue="Fragment",
    binwidth=2,
    alpha=0.5,
    palette=palette
)

plt.xticks(range(0, int(barcode_mut["barcode_per_mut"].max()) + 10, 10))
plt.xlabel("Number of barcode per mutation (AA)")
plt.ylabel("Count")
plt.show()

## PLOT nb barcode per AA
barcode_mut = barcode_mut[(barcode_mut['Fragment'] == 'F13') | (barcode_mut['Fragment'] == 'F43')]
fragments = barcode_mut['Fragment'].unique()
num_fragments = len(fragments)
nplots = 1

for start in range(0, num_fragments, nplots):
    subset = fragments[start:start+nplots]
    
    # Define grid size
    cols = math.ceil(math.sqrt(len(subset)))
    rows = math.ceil(len(subset) / cols)

    # Create subplots
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3.9), constrained_layout=True, dpi=1000)
    axes = axes.flatten() if len(subset) > 1 else [axes]
    
    for i, fragment in enumerate(subset):
        ax = axes[i]
        ax.set_aspect(1)
    
        # Subset data for the fragment
        df_subset = barcode_mut[barcode_mut['Fragment'] == fragment].pivot(
            index='mutation_aa', columns='position', values='barcode_per_mut'
        )
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
            cbar=True, vmin = vmin, vmax = vmax, norm = norm
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
        
            # Check if the WT aa (ref_codon) exists in the Y-axis (mutations)
            if ref_aa in df_subset.index and mut_aa in df_subset.columns:
                y_pos = df_subset.index.get_loc(ref_aa)  # Get row index
                x_pos = df_subset.columns.get_loc(mut_aa)  # Get column index
                
                # ✅ Add gray square for WT codons
                ax.add_patch(plt.Rectangle((x_pos + 0.05, y_pos + 0.05), 0.9, 0.9, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Show plots
    plt.show()
    # first_fragment = subset[0]
    # last_fragment = subset[-1]
    # last_fragment_num = last_fragment.replace('F','')
    # plt.savefig(f"{wkdir}/figures_reviews/golden_gate_nb_barcodes_{first_fragment}.png")

#%% Compute %barcode diversity AA
df_ref['expected_mutations'] = df_ref['Fragment'].map(
    lambda frag: 20 * (frag_len[frag] / 3)  # 19 possible amino acid changes per position + stop 
)

# Merge updated reference table with barcode_mut
barcode_mut = pd.merge(barcode_mut, df_ref[['Fragment', 'expected_mutations']]).drop_duplicates()
barcode_mut = barcode_mut[barcode_mut['mutation_aa'] != barcode_mut['ref_aa']]

# Create binary columns based on the threshold values
barcode_mut['barcode_0'] = (barcode_mut['barcode_per_mut'] > 0).astype(int)
barcode_mut['barcode_4'] = (barcode_mut['barcode_per_mut'] > 4).astype(int)
barcode_mut['barcode_9'] = (barcode_mut['barcode_per_mut'] > 9).astype(int)

# Group by Fragment, aggregate the sums and counts
dg = barcode_mut.groupby('Fragment').agg(
    At_least_1=('barcode_0', 'sum'),
    More_than_4=('barcode_4', 'sum'),
    More_than_9=('barcode_9', 'sum'),
    n=('expected_mutations', 'first')
).reset_index()

# Calculate the percentages
dg['>0'] = (dg['At_least_1'] / dg['n']) * 100
dg['>4'] = (dg['More_than_4'] / dg['n']) * 100
dg['>9'] = (dg['More_than_9'] / dg['n']) * 100

df = dg[['Fragment','>0','>4','>9']]
df.to_csv('/home/alicia-pageau/Documents/antifungal_project/PDR1/golden_gate/05_analyze/golden_gate_barcode_diversity.csv')
