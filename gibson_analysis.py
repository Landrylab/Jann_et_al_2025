#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:02:31 2025

@author: alicia

17-02-2025
PDR1 Gibson Sequencing Analysis
"""

import pandas as pd
import pyarrow
import numpy as np
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
import glob
import os
import re
import itertools
from itertools import islice
import math
from natsort import natsorted
import matplotlib.patches as mpatches

print('pandas', pd.__version__)
print('matplotlib', plt.matplotlib.__version__)
print('numpy', np.__version__)
print('seaborn', sns.__version__)

wkdir = '/home/alicia/Documents/antifungal_project/PDR1'
os.chdir(wkdir)

# aa order per properties for heatmaps
aa_order = ["*", "P", "G", "C", "Q", "N", "T", "S", "E", "D",
            "K", "H", "R", "W", "Y", "F", "M", "L", "I", "V", "A"]


#%% Functions
def first_mutation_position(sequence1, sequence2):
    pos = [i for i,x in enumerate(zip(sequence1,sequence2)) if x[0]!=x[1]]
    return pos[0]

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

def is_valid_nnk(codon):
    return bool(re.match(r'^[ACGT]{2}[GT]$', codon))  # First two bases: A/C/G/T, Last base: G/T

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

R = {rec.id.split('_')[1]: str(rec.seq).upper() for rec in dref}
D = {rec.id.split('_')[1]:len(rec.seq) for rec in dref}

# Make a df with reference codon (WT)
codon_data = []
for seq_id, sequence in R.items():
    for i in range(0, len(sequence) - 2, 3):  # Iterate in steps of 3
        codon = sequence[i:i+3]  # Extract codon
        codon_position = i // 3 + 1  # Convert base position to codon number
        codon_data.append([seq_id, codon_position, codon])
df_ref = pd.DataFrame(codon_data, columns=["Fragment", "mut_codon", "ref_codon"])
df_ref['position'] = ((df_ref['Fragment'].str.replace("F","").astype(int) -1)*25) + df_ref['mut_codon']

#%% Processing reads
reads = glob.glob(f"{wkdir}/gibson/04_merged/01_2025/Gibson*.fasta")

STATS = {}
M = []

for fread in reads:
    # Read fasta file
    print(fread)
    fragment = fread.split('/')[-1].split('.')[0].replace("Gibson","")
    fragment_len = tot_len[fragment]
    dread = [str(rec.seq) for rec in SeqIO.parse(fread,"fasta")]

    # Reading and processing of reads
    dR = pd.DataFrame({"seq" : dread})
    a = dR.shape[0]
    print(dR.shape)
    if dR.shape[0] > 0:
        # Removing unaligned
        dR['len'] = dR['seq'].apply(lambda x: len(x))
        dR = dR[dR['len'] == fragment_len].reset_index(drop=True)
        b = dR.shape[0]
        print(dR.shape)
        if dR.shape[0] > 0:
            # Subset core and barcode
            dR['core'] = dR['seq'].apply(lambda x: x[extra_front[fragment]:-barcode_length-extra_back[fragment]])
            dR['barcode'] = dR['seq'].apply(lambda x: x[-barcode_length:])
            dR['pattern'] = dR['barcode'].apply(lambda x: x[2:14] + "_" + x[16:28])
            # Calculating mismatches
            dR['ref'] = R[fragment]
            dR['mismatches'] = dR.apply(lambda x: len([(a,b) for a,b in zip(x['core'],x['ref']) if a!=b]), axis=1)
            
            # Removing reads with more than 1 codon mutated
            dR['ref_codon'] = dR.apply(lambda x: check_codons(x['core'], x['ref'])[0], axis=1)
            dR['mutation'] = dR.apply(lambda x: check_codons(x['core'], x['ref'])[1], axis=1)
            dR['changed_codons_n'] = dR.apply(lambda x: check_codons(x['core'], x['ref'])[2], axis=1)
            # Include single mutation and a wildtype
            dR = dR[dR['changed_codons_n'] <= 1].reset_index(drop=True)
            c = dR.shape[0]
            print(dR.shape)
            if dR.shape[0] > 0:
            
                # Removing sequences with N in the sequence
                dR = dR[~dR['core'].str.contains('N')].reset_index(drop=True)
                dR = dR[~dR['barcode'].str.contains('N')].reset_index(drop=True)
                print(dR.shape)
                d = dR.shape[0]
                dR.head()

                if dR.shape[0] > 0:
                    # Get mutated positions
                    dR['mut_pos'] = dR.apply(lambda x: mutation_positions(x['core'],x['ref'])[0], axis=1)
                    dR['mut_codon'] = dR['mut_pos'].apply(lambda x: int(np.floor((x-1)/3))+1 if x>0 else 0)
                    
                    # Merging reads with barcode patterns
                    dM = pd.merge(dR, dbarc, on = ['pattern'], how = 'left')
                    # Some patterns are not matching (because of mismatches)
                    nonmatching_barcodes = dM[dM['codon'].isnull()].shape[0]
                    print("Number of nonmatching barcodes %d" %(nonmatching_barcodes))

                    dM = dM[~dM['codon'].isnull()].reset_index(drop=True)
                    dM['Fragment'] = fragment
                    
                    dM.to_csv(f"{wkdir}/gibson/05_analyze/01_2025/"+fragment+".csv", sep=",", header=True, index=False)
                    M.append(dM)
                    
                    e = dM.shape[0]
                    print(dM.shape)
                    STATS[fragment] = [a,b,c,d,e]

dSTAT = pd.DataFrame(STATS).T.reset_index()
dSTAT.columns = ['fragment','N_reads','N_aligned','N_1codon','N_noN','N_barcode']
dSTAT['retained'] = dSTAT.apply(lambda x: (x['N_barcode']/x['N_reads'])*100, axis=1)
dSTAT = dSTAT.sort_values(by=['fragment']).reset_index(drop=True)
dSTAT.to_csv(f"{wkdir}/gibson/05_analyze/01_2025/gibson_STATS.tsv", sep='\t', header=True, index=False)

#%% Stats per fragment
read_depth_threshold = 1

result_files = glob.glob(f"{wkdir}/gibson/05_analyze/01_2025/F*.csv") + \
               glob.glob(f"{wkdir}/gibson/05_analyze/F*.csv")

# Define regex patterns for expected groups
expected_patterns = [
    [r"F1_.*\.csv", r"F[2-9]\.csv"],   
    [r"F1[0-8]\.csv", r"F13_.*\.csv"],  
    ["F19.csv",r"F2[0-7]\.csv"],
    [r"F2[8-9]\.csv", r"F3[0-6]\.csv"],  
    [r"F3[7-9]\.csv", r"F4[0-2]\.csv", r"F43_.*\.csv"]  
]

all_results = []
all_informative_barcode_stat = []

# Process files that match regex pattern
for pattern_group in expected_patterns:
    # Filter files that match ANY pattern in the group
    batch_files = [f for f in result_files if any(re.search(p, os.path.basename(f)) for p in pattern_group)]
    batch_files = natsorted(batch_files, key=os.path.basename)
    rM = pd.concat((pd.read_csv(f) for f in batch_files), ignore_index=True)
    all_results.append(rM)
    
    read_depth = rM.groupby(['Fragment','barcode','core']).agg(read_depth = ('seq','count')).reset_index()
    read_depth = read_depth[read_depth['read_depth'] > read_depth_threshold]

    rM_read_depth = pd.merge(rM, read_depth, on = ['Fragment','barcode','core'])

    ## Barcode uniqueness - Get barcode associated to only one mutation
    # Select unique barcode
    barcode_count = rM_read_depth.groupby(['Fragment', 'barcode']).agg(barcode_reads = ('core','nunique')).reset_index()
    unique_barcode = barcode_count[barcode_count['barcode_reads'] == 1].drop(columns='barcode_reads')
    rM_unique_barcode = pd.merge(rM_read_depth, unique_barcode, on = ['Fragment','barcode'])
    # Calculate % informative barcodes
    pourcent_informative_barcode = (
        barcode_count[barcode_count['barcode_reads'] == 1]
        .groupby('Fragment')
        .size()  # Count barcodes with unique core
        / barcode_count.groupby('Fragment').size()  # Total barcodes per Fragment
    ) * 100
    pourcent_informative_barcode = pourcent_informative_barcode.reset_index(name='%_informative_barcode')
    
    # Count informative barcodes (all)
    nb_informative_barcode = (
        barcode_count[barcode_count['barcode_reads'] == 1]
        .groupby('Fragment')
        .agg(nb_informative_barcode=('barcode_reads', 'sum')).reset_index()
    )
    
    # Calculate % informative barcodes excluding WT (mutation == 'REF')
    barcode_count_mut = pd.merge(rM_read_depth[['Fragment','barcode','mutation']].drop_duplicates(), barcode_count, on = ['Fragment','barcode'])    
    pourcent_informative_barcode_excl_WT = (
        barcode_count_mut[(barcode_count_mut['barcode_reads'] == 1) & (barcode_count_mut['mutation'] != 'REF')]
        .groupby('Fragment')
        .size()  # Count barcodes with unique core, excluding WT
        / barcode_count.groupby('Fragment').size()  # Total barcodes per Fragment excluding WT
    ) * 100
    
    pourcent_informative_barcode_excl_WT = pourcent_informative_barcode_excl_WT.reset_index(name='%_informative_barcode_excl_WT')
    
    # Count informative barcodes excluding WT (mutation == 'REF')
    nb_informative_barcode_excl_WT = (
        barcode_count_mut[(barcode_count_mut['barcode_reads'] == 1) & (barcode_count_mut['mutation'] != 'REF')]
        .groupby('Fragment')
        .agg(nb_informative_barcode_excl_WT=('barcode_reads', 'sum')).reset_index()
    )

    ## Get % WT reads
    wt_count = rM_read_depth.groupby(['Fragment','mutation']).agg(read_per_mutation = ('mutation','count')).reset_index()
    pourcent_wt = (
        wt_count[wt_count['mutation'] == 'REF']
        .groupby('Fragment')['read_per_mutation']
        .sum()  # Sum WT reads per fragment
        / wt_count.groupby('Fragment')['read_per_mutation'].sum()  # Total reads per fragment
    ) * 100
    pourcent_wt = pourcent_wt.reset_index(name='WT_read_percentage')
    
    # Merge % informative barcode with % informative barcode excluding WT
    informative_barcode_stat = pourcent_informative_barcode.merge(
        nb_informative_barcode, on='Fragment'
    ).merge(
        pourcent_informative_barcode_excl_WT, on='Fragment'
    ).merge(
        nb_informative_barcode_excl_WT, on='Fragment'
    ).merge(
        pourcent_wt, on='Fragment'
    )
    all_informative_barcode_stat.append(informative_barcode_stat)

stats = pd.concat(all_informative_barcode_stat, ignore_index=True)
stats['%_barcode_WT'] = stats['%_informative_barcode'] - stats['%_informative_barcode_excl_WT']

# Plot
sns.set_theme(style="whitegrid")
plt.figure(figsize=(12, 8), dpi=300)
sns.barplot(x="Fragment", y="%_informative_barcode", data=stats, color='lightblue')
sns.barplot(x="Fragment", y="%_barcode_WT", data=stats, color='red')
plt.xticks(rotation=90)
top_bar = mpatches.Patch(color='lightblue', label='Mutations')
bottom_bar = mpatches.Patch(color='red', label='Wild-type')
plt.legend(handles=[top_bar, bottom_bar],loc="upper right")
plt.show()


#%% Figures (codon)
read_depth_threshold = 1

result_files = glob.glob(f"{wkdir}/gibson/05_analyze/F43_4.csv") +\
               glob.glob(f"{wkdir}/gibson/05_analyze/F43_16.csv") +\
               glob.glob(f"{wkdir}/gibson/05_analyze/F13_[3-5].csv") +\
               glob.glob(f"{wkdir}/gibson/05_analyze/F13_13.csv")
rM = pd.concat((pd.read_csv(f) for f in result_files), ignore_index=True)

read_depth = rM.groupby(['Fragment','barcode','core']).agg(read_depth = ('seq','count')).reset_index()
read_depth = read_depth[read_depth['read_depth'] > read_depth_threshold]

rM_read_depth = pd.merge(rM, read_depth, on = ['Fragment','barcode','core'])

barcode_count = rM_read_depth.groupby(['Fragment', 'barcode']).agg(barcode_reads = ('core','nunique')).reset_index()
unique_barcode = barcode_count[barcode_count['barcode_reads'] == 1].drop(columns='barcode_reads')
rM_unique_barcode = pd.merge(rM_read_depth, unique_barcode, on = ['Fragment','barcode'])

# Barcode diversity - Nb unique barcode per mutation
barcode_mut = rM_unique_barcode.groupby(['Fragment','core','mutation','mut_codon','ref_codon']).agg(barcode_per_mut = ('barcode','nunique'), reads = ('seq','count')).reset_index()
barcode_mut = barcode_mut[barcode_mut['mutation'].apply(is_valid_nnk)]

# Add position mutation in the whole protein seq
barcode_mut['position'] = ((barcode_mut['Fragment'].str.replace("F","").astype(int) -1)*25) + barcode_mut['mut_codon']

# Make a df with reference codon (WT)
codon_data = []
for seq_id, sequence in R.items():
    for i in range(0, len(sequence) - 2, 3):  # Iterate in steps of 3
        codon = sequence[i:i+3]  # Extract codon
        codon_position = i // 3 + 1  # Convert base position to codon number
        codon_data.append([seq_id, codon_position, codon])
df_ref = pd.DataFrame(codon_data, columns=["Fragment", "mut_codon", "ref_codon"])
df_ref['position'] = ((df_ref['Fragment'].str.replace("F","").astype(int) -1)*25) + df_ref['mut_codon']

# Create a full DataFrame with all combinations
fragments = set(barcode_mut['Fragment'].unique())  # Convert to set for faster lookup
full_index = [(frag, codon) for frag, length in frag_len.items() if frag in fragments for codon in range(1, (length//3) + 1)]
full_df = pd.DataFrame(full_index, columns=['Fragment', 'mut_codon'])
full_df = pd.merge(full_df, df_ref, how = 'left')

barcode_mut = full_df.merge(barcode_mut, how='left')
barcode_mut['Fragment'] = pd.Categorical(barcode_mut['Fragment'], categories=natsorted(barcode_mut['Fragment'].unique()), ordered=True)
barcode_mut = barcode_mut.sort_values(['Fragment', 'mut_codon'])

#### PLOT nb barcode per mutation
# Get unique fragments
fragments = barcode_mut['Fragment'].unique()
num_fragments = len(fragments)

# Define grid size
cols = math.ceil(math.sqrt(num_fragments))  
rows = math.ceil(num_fragments / cols)

# Create subplots
fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5), constrained_layout=True)
axes = axes.flatten() if num_fragments > 1 else [axes]

# Set color scale limits
vmin, vmax = 0, 10  

for i, fragment in enumerate(fragments):
    ax = axes[i]
    ax.set_aspect(1)

    # Subset data for the fragment
    df_subset = barcode_mut[barcode_mut['Fragment'] == fragment].pivot(
        index='mutation', columns='position', values='barcode_per_mut'
    )
    
    # Plot heatmap with fixed color scale
    sns.heatmap(
        df_subset, cmap="Blues", annot=False, linewidths=0.5, ax=ax,
        vmin=vmin, vmax=vmax, cbar=True
    )

    # Titles and labels
    ax.set_title(f"{fragment}")
    ax.set_xlabel("Mutated Codon")
    ax.set_ylabel("Mutation")

    # ✅ Adjust axis labels
    ax.set_yticks(np.arange(len(df_subset.index)) + 0.5)
    #ax.set_yticklabels(df_subset.index, rotation=0)

    xtick_positions = list(range(0, len(df_subset.columns), 14))
    last_pos = int(len(df_subset.columns) - 1)
    if last_pos not in xtick_positions:
        xtick_positions.append(last_pos)  # Append correctly
    xtick_positions = np.array(xtick_positions, dtype=int)  # Ensure integer array
    xtick_labels = df_subset.columns[xtick_positions].tolist()  # Ensure a list format
    ax.set_xticks(xtick_positions + 0.5)  # Centered tick positions
    ax.set_xticklabels(xtick_labels, rotation=0)  # Assign labels
    
    # ✅ Add gray squares for WT codons
    for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
        ref_codon = str(row['ref_codon'])  # Ensure it's a string
        mut_codon = row['position']  # This is a position (number)
    
        # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
        if ref_codon in df_subset.index and mut_codon in df_subset.columns:
            y_pos = df_subset.index.get_loc(ref_codon)  # Get row index
            x_pos = df_subset.columns.get_loc(mut_codon)  # Get column index
            
            # ✅ Add gray square for WT codons
            ax.add_patch(plt.Rectangle((x_pos + 0.05, y_pos +0.05), 0.95, 0.95, color='gray', lw=0.1))

# Remove unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Show plots
plt.show()

#### PLOT %barcode diversity
from itertools import product
N = ['A', 'T', 'C', 'G']
K = ['G', 'T']

# Generate all NNK codons
nnk_codons = {''.join(codon) for codon in product(N, N, K)}

# Convert to dictionary for easy lookup
nnk_dict = {codon: None for codon in nnk_codons}

# Function to count WT-encoding codons per position
def count_wt_nnk(ref_codon):
    # Find how many NNK codons encode the same amino acid as ref_codon
    return sum(1 for codon in nnk_dict if codon == ref_codon)

# Apply WT-encoding count per fragment
df_ref['wt_nnk_count'] = df_ref['ref_codon'].apply(count_wt_nnk)

# Compute total WT-encoding NNK mutations per fragment
wt_nnk_counts = df_ref.groupby('Fragment')['wt_nnk_count'].sum().reset_index()
wt_nnk_counts = wt_nnk_counts.rename(columns={'wt_nnk_count':'wt_nnk_count_per_frag'})
df_ref2 = pd.merge(df_ref,wt_nnk_counts)

# Compute expected number of mutations per fragment
df_ref2['expected_mutations'] = df_ref2['Fragment'].map(
    lambda frag: 32 * (frag_len[frag]/3)  # 32 possible NNK codons per position
)

# Compute total expected mutations per fragment
df_ref2['nb_possible_mutation'] = df_ref2['expected_mutations'] - df_ref2['wt_nnk_count_per_frag']

barcode_mut = pd.merge(barcode_mut, df_ref2[['Fragment','nb_possible_mutation']]).drop_duplicates()
# Create binary columns based on the threshold values
barcode_mut['barcode_0'] = (barcode_mut['barcode_per_mut'] > 0).astype(int)
barcode_mut['barcode_4'] = (barcode_mut['barcode_per_mut'] > 4).astype(int)
barcode_mut['barcode_9'] = (barcode_mut['barcode_per_mut'] > 9).astype(int)

# Group by Fragment, aggregate the sums and counts
dg = barcode_mut.groupby('Fragment').agg(
    At_least_1=('barcode_0', 'sum'),
    More_than_4=('barcode_4', 'sum'),
    More_than_9=('barcode_9', 'sum'),
    n=('nb_possible_mutation', 'first')
).reset_index()

# Calculate the percentages
dg['>0'] = (dg['At_least_1'] / dg['n']) * 100
dg['>4'] = (dg['More_than_4'] / dg['n']) * 100
dg['>9'] = (dg['More_than_9'] / dg['n']) * 100


df = dg[['Fragment','>0','>4','>9']]

# Set 'Fragment' as the column index for proper heatmap structure
df_melted = df.melt(id_vars=['Fragment'], var_name='nb_barcode', value_name='Percentage')
df_melted['Fragment'] = pd.Categorical(df_melted['Fragment'], categories=natsorted(df_melted['Fragment'].unique()), ordered=True)
df_melted = df_melted.sort_values(['Fragment'])

# Pivot table to reshape for heatmap (Threshold as rows, Fragment as columns)
df_pivot = df_melted.pivot(index='nb_barcode', columns='Fragment', values='Percentage')

# 🔥 Create heatmap
plt.figure(figsize=(10, 4))
ax = sns.heatmap(df_pivot, annot=True, cmap="Blues", fmt=".1f", vmin=0, vmax=100)
ax.set_aspect(1)
# 🎯 Formatting
ax.set_xlabel("Fragment")
ax.set_ylabel("")
ax.set_title("% barcode diversity")

plt.show()

#### PLOT read_depth per mutation
# Get unique fragments
fragments = barcode_mut['Fragment'].unique()
num_fragments = len(fragments)

# Define grid size
cols = math.ceil(math.sqrt(num_fragments))  
rows = math.ceil(num_fragments / cols)

# Create subplots
fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5), constrained_layout=True)
axes = axes.flatten() if num_fragments > 1 else [axes]

# Set color scale limits
vmin, vmax = 0, 2000

for i, fragment in enumerate(fragments):
    ax = axes[i]

    # Subset data for the fragment
    df_subset = barcode_mut[barcode_mut['Fragment'] == fragment].pivot(
        index='mutation', columns='position', values='reads'
    )

    # Plot heatmap with fixed color scale
    sns.heatmap(
        df_subset, cmap="Reds", annot=False, linewidths=0.5, ax=ax,
        vmin=vmin, vmax=vmax, cbar=True
    )

    # Titles and labels
    ax.set_title(f"{fragment}")
    ax.set_xlabel("Mutated Codon")
    ax.set_ylabel("Mutation")

    # ✅ Adjust axis labels
    ax.set_yticks(np.arange(len(df_subset.index)) + 0.5)
    #ax.set_yticklabels(df_subset.index, rotation=0)
    
    # ✅ Add gray squares for WT codons
    for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
        ref_codon = str(row['ref_codon'])  # Ensure it's a string
        mut_codon = row['position']  # This is a position (number)
    
        # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
        if ref_codon in df_subset.index and mut_codon in df_subset.columns:
            y_pos = df_subset.index.get_loc(ref_codon)  # Get row index
            x_pos = df_subset.columns.get_loc(mut_codon)  # Get column index
            
            # ✅ Add gray square for WT codons
            ax.add_patch(plt.Rectangle((x_pos + 0.05, y_pos +0.05), 0.95, 0.95, color='gray', lw=0.1))

# Remove unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Show plots
plt.show()

#### PLOT coverage per codon
# Compute coverage of nucleic acide mutation at every position --> should be 3,3,1 for NNK
#dcov_mut = barcode_mut.groupby(['changed_codons_n', 'mut_pos','Fragment']).agg(n_muts = ('nt','nunique')).reset_index()

# Compute coverage of mutation at every codon --> should be 32 for NNK
dcov_codon = barcode_mut.groupby(['mut_codon', 'Fragment']).agg(n_muts = ('mutation','count')).reset_index()

plt.figure(figsize=(10, 4))

heatmap_data = dcov_codon.pivot(index="Fragment", columns="mut_codon", values="n_muts")
sns.heatmap(heatmap_data, cmap="Blues", annot=True, linewidths=0.5,
    vmin=0, vmax=32, cbar=True)

ax.set_title("Codon coverage")
ax.set_xlabel("Mutated Codon")
ax.set_ylabel("Fragments")

plt.show()


#%% Figures (AA)
read_depth_threshold = 1

# Get result files
result_files = glob.glob(f"{wkdir}/gibson/05_analyze/F43_4.csv") +\
               glob.glob(f"{wkdir}/gibson/05_analyze/F43_16.csv") +\
               glob.glob(f"{wkdir}/gibson/05_analyze/F13_[3-5].csv") +\
               glob.glob(f"{wkdir}/gibson/05_analyze/F13_13.csv")


rM = pd.concat((pd.read_csv(f) for f in result_files), ignore_index=True)

read_depth = rM.groupby(['Fragment','barcode','core']).agg(read_depth = ('seq','count')).reset_index()
read_depth = read_depth[read_depth['read_depth'] > read_depth_threshold]

rM_read_depth = pd.merge(rM, read_depth, on = ['Fragment','barcode','core'])

barcode_count = rM_read_depth.groupby(['Fragment', 'barcode']).agg(barcode_reads = ('core','nunique')).reset_index()
unique_barcode = barcode_count[barcode_count['barcode_reads'] == 1].drop(columns='barcode_reads')
rM_unique_barcode = pd.merge(rM_read_depth, unique_barcode, on = ['Fragment','barcode'])

# Barcode diversity - Nb unique barcode per AA

rM_unique_barcode['mutation_aa'] = rM_unique_barcode['mutation'].map(codon_to_aa)
rM_unique_barcode['ref_aa'] = rM_unique_barcode['ref_codon'].map(codon_to_aa)

barcode_mut = rM_unique_barcode.groupby(['Fragment','mutation_aa','mut_codon']).agg(barcode_per_mut = ('barcode','nunique')).reset_index()
barcode_mut = barcode_mut[barcode_mut['mutation_aa'] != 'REF']


# Add position mutation in the whole protein seq
barcode_mut['position'] = ((barcode_mut['Fragment'].str.replace("F","").astype(int) -1)*25) + barcode_mut['mut_codon']

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

# Create a full DataFrame with all combinations
fragments = set(barcode_mut['Fragment'].unique())  # Convert to set for faster lookup
full_index = [(frag, codon) for frag, length in frag_len.items() if frag in fragments for codon in range(1, (length//3) + 1)]
full_df = pd.DataFrame(full_index, columns=['Fragment', 'mut_codon'])
full_df = pd.merge(full_df, df_ref, how = 'left')

barcode_mut = full_df.merge(barcode_mut, how='left')
barcode_mut['Fragment'] = pd.Categorical(barcode_mut['Fragment'], categories=natsorted(barcode_mut['Fragment'].unique()), ordered=True)
barcode_mut = barcode_mut.sort_values(['Fragment', 'mut_codon'])

# PLOT nb barcode per mutation
# Get unique fragments
fragments = barcode_mut['Fragment'].unique()
num_fragments = len(fragments)

# Define grid size
cols = math.ceil(math.sqrt(num_fragments))  
rows = math.ceil(num_fragments / cols)

# Create subplots
fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3.9), constrained_layout=True, dpi=1000)
axes = axes.flatten() if num_fragments > 1 else [axes]

# Set color scale limits
vmin, vmax = 0, 10

for i, fragment in enumerate(fragments):
    ax = axes[i]
    ax.set_aspect(1)

    # Subset data for the fragment
    df_subset = barcode_mut[barcode_mut['Fragment'] == fragment].pivot(
        index='mutation_aa', columns='position', values='barcode_per_mut'
    )
    df_subset = df_subset.loc[aa_order]

    # Plot heatmap with fixed color scale
    sns.heatmap(
        df_subset, cmap="Blues", annot=False, linewidths=0.5, ax=ax,
        vmin=vmin, vmax=vmax, cbar=True
    )

    # Titles and labels
    ax.set_title(f"{fragment}")
    ax.set_xlabel("Codon position")
    ax.set_ylabel("")

    # ✅ Adjust axis labels
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
    
    # ✅ Add gray squares for WT codons
    for idx, row in df_ref[df_ref['Fragment'] == fragment].iterrows():
        ref_codon = str(row['ref_aa'])  # Ensure it's a string
        mut_codon = row['position']  # This is a position (number)
    
        # Check if the WT codon (ref_codon) exists in the Y-axis (mutations)
        if ref_codon in df_subset.index and mut_codon in df_subset.columns:
            y_pos = df_subset.index.get_loc(ref_codon)  # Get row index
            x_pos = df_subset.columns.get_loc(mut_codon)  # Get column index
            
            # ✅ Add gray square for WT codons
            ax.add_patch(plt.Rectangle((x_pos + 0.05, y_pos +0.05), 0.95, 0.95, color='gray', lw=0.1))

# Remove unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Show plots
plt.show()


# PLOT %barcode diversity AA

df_ref['expected_mutations'] = df_ref['Fragment'].map(
    lambda frag: 20 * (frag_len[frag] / 3)  # 19 possible amino acid changes per position
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

# Set 'Fragment' as the column index for proper heatmap structure
df_melted = df.melt(id_vars=['Fragment'], var_name='nb_barcode', value_name='Percentage')
df_melted['Fragment'] = pd.Categorical(df_melted['Fragment'], categories=natsorted(df_melted['Fragment'].unique()), ordered=True)
df_melted = df_melted.sort_values(['Fragment'])

# Pivot table to reshape for heatmap (Threshold as rows, Fragment as columns)
df_pivot = df_melted.pivot(index='nb_barcode', columns='Fragment', values='Percentage')

# 🔥 Create heatmap
plt.figure(figsize=(10, 4), dpi=300)
ax = sns.heatmap(df_pivot, annot=True, cmap="Blues", fmt=".1f", vmin=0, vmax=100)
ax.set_aspect(1)

# 🎯 Formatting
ax.set_xlabel("Fragment")
ax.set_ylabel("")
ax.set_title("Barcode diversity (%)")

plt.show()
