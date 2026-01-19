#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:02:31 2025

@author: alicia-pageau

10-11-2025
PDR1 Gibson Sequencing Analysis
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
dbarc['mut_codon'] = dbarc['codon'] - dbarc.groupby('Fragment')['codon'].transform('first') + 1

dref = list(SeqIO.parse("references_PDR1_corrected.fasta","fasta"))

R = {rec.id.split('_')[1]: str(rec.seq).upper() for rec in dref}
D = {rec.id.split('_')[1]:len(rec.seq) for rec in dref}

# Make a df with reference (WT)
codon_data = []
for seq_id, sequence in R.items():
    for i in range(0, len(sequence) - 2, 3):  # Iterate in steps of 3
        codon = sequence[i:i+3]  # Extract codon
        codon_position = i // 3 + 1  # Convert base position to codon number
        amino_acid = codon_to_aa.get(codon) # Convert codon to amino acid 
        codon_data.append([seq_id, codon_position, codon, amino_acid])
df_ref = pd.DataFrame(codon_data, columns=["Fragment", "mut_codon", "ref_codon", "ref_aa"])
df_ref['position'] = ((df_ref['Fragment'].str.replace("F","").astype(int) -1)*25) + df_ref['mut_codon']
df_ref['expected_mutations'] = df_ref['Fragment'].map(
    lambda frag: 32 * (frag_len[frag]/3)  # 32 possible NNK codons per position
)
#%% Processing reads
reads = glob.glob(f"{wkdir}/gibson/04_merged/01_2025/GibsonF*.fasta.gz") + \
    glob.glob(f"{wkdir}//gibson/04_merged/10_2025/Gibson_F38_V2.fasta.gz") + \
    glob.glob(f"{wkdir}/gibson/04_merged/PDR1_GiF1_*") + \
    glob.glob(f"{wkdir}/gibson/04_merged/PDR1_GiF13_[3-5]*") + \
    glob.glob(f"{wkdir}/gibson/04_merged/PDR1_GiF13_13*") + \
    glob.glob(f"{wkdir}/gibson/04_merged/PDR1_GiF43_4*") + \
    glob.glob(f"{wkdir}/gibson/04_merged/PDR1_GiF43_16*")

STATS = {}
M = []

for fread in reads:
    # Read fasta file
    print(fread)
    fragment = fread.split('/')[-1].split('.')[0].replace("Gibson","")
    # fragment = fread.split('/')[-1].split('.')[0].split('_')[1]
    # fragment = fread.split('/')[-1].split('.')[0].split('_')[1].replace("Gi","")
    # pcr = fread.split('/')[-1].split('.')[0].split('_')[-1]

    print(fragment)
    fragment_len = tot_len[fragment]
    
    open_func = gzip.open if fread.endswith(".gz") else open
    with open_func(fread, "rt") as handle:
       dread = [str(rec.seq) for rec in SeqIO.parse(handle, "fasta")]
    

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
                    dM = pd.merge(dR, dbarc, on='pattern', how = 'left')
                    # Some patterns are not matching (because of mismatches)
                    nonmatching_barcodes = dM[dM['codon'].isnull()].shape[0]
                    print("Number of nonmatching barcodes %d" %(nonmatching_barcodes))
                    dM = dM[~dM['codon'].isnull()].reset_index(drop=True)
                    dM['Fragment'] = fragment
                    
                    dM.to_csv(f"{wkdir}/gibson/05_analyze/"+fragment+".csv", sep=",", header=True, index=False)
                    M.append(dM)
                    
                    e = dM.shape[0]
                    print(dM.shape)
                    
                    STATS[fragment] = [a,b,c,d,e]

dSTAT = pd.DataFrame(STATS).T.reset_index()
dSTAT.columns = ['fragment','N_reads','N_aligned','N_1codon','N_noN','N_barcode']
dSTAT['retained'] = dSTAT.apply(lambda x: (x['N_barcode']/x['N_reads'])*100, axis=1)
dSTAT = dSTAT.sort_values(by=['fragment']).reset_index(drop=True)
dSTAT.to_csv(f"{wkdir}/gibson/05_analyze/gibson_STATS.tsv", sep='\t', header=True, index=False)

#%% Combine results
read_depth_threshold = 1
gibson_result_files = glob.glob(f"{wkdir}/gibson/05_analyze/F1_*.csv.gz") + \
               glob.glob(f"{wkdir}/gibson/05_analyze/F13_[3-5].csv.gz") + \
               glob.glob(f"{wkdir}/gibson/05_analyze/F13_13.csv.gz") + \
               glob.glob(f"{wkdir}/gibson/05_analyze/F43_4.csv.gz") + \
               glob.glob(f"{wkdir}/gibson/05_analyze/F43_16.csv.gz") + \
               glob.glob(f"{wkdir}/gibson/05_analyze/01_2025/*.csv.gz")  + \
               glob.glob(f"{wkdir}/gibson/05_analyze/10_2025/F38.csv.gz")
gibson_result_files = [f for f in gibson_result_files if "F38_v1" not in f]
dfs = []
for f in gibson_result_files:
    df = pd.read_csv(f)
    df['read_depth'] = df.groupby(['Fragment', 'barcode', 'core'])['seq'].transform('count')
    df = df[df['read_depth'] > read_depth_threshold].drop_duplicates()
    dfs.append(df)
gibson = pd.concat(dfs, ignore_index=True)
#gibson.to_csv(f"{wkdir}/gibson_all_fragments_sequencing_read_2010205.csv", index=False)
#%% Stats per fragment
rM_read_depth = pd.read_csv(f"{wkdir}/gibson_all_fragments_sequencing_read_20102025.csv")

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

informative_barcode_stat['%_barcode_WT'] = informative_barcode_stat['%_informative_barcode'] - informative_barcode_stat['%_informative_barcode_excl_WT']
informative_barcode_stat.to_csv(f"{wkdir}/gibson/05_analyze/gibson_informative_barcodes_stats.csv")


#%% Compute uniformity (codon-level)
rM_read_depth=pd.read_csv(f"{wkdir}/gibson_all_fragments_sequencing_read_20102025.csv")

barcode_count = rM_read_depth.groupby(['Fragment', 'barcode']).agg(barcode_reads = ('core','nunique')).reset_index()
unique_barcode = barcode_count[barcode_count['barcode_reads'] == 1].drop(columns='barcode_reads')
rM_unique_barcode = pd.merge(rM_read_depth, unique_barcode, on = ['Fragment','barcode'])

# Barcode diversity - Nb unique barcode per mutation
barcode_mut = rM_unique_barcode.groupby(['Fragment','core','mutation','mut_codon','ref_codon']).agg(barcode_per_mut = ('barcode','nunique'), reads = ('seq','count')).reset_index()
barcode_mut = barcode_mut[barcode_mut['mutation'].apply(is_valid_nnk)]
barcode_mut['read_freq'] = barcode_mut['reads'] / barcode_mut.groupby('Fragment')['reads'].transform('sum')
barcode_mut['mutation_aa'] = barcode_mut['mutation'].map(codon_to_aa)
barcode_mut['ref_aa'] = barcode_mut['ref_codon'].map(codon_to_aa)
barcode_mut = barcode_mut.fillna(0)

# # In prevision of plot with aa as well
# aa_counts = (
#     barcode_mut
#     .groupby(['position', 'mutation_aa'])
#     .agg(barcode_per_mut_aa=('barcode_per_mut', 'sum'))
#     .reset_index()
# )
# barcode_mut = barcode_mut.merge(aa_counts, on=['position', 'mutation_aa'], how='left')
# barcode_mut.to_csv('/home/alicia-pageau/Documents/antifungal_project/PDR1/00_scripts/Jann_et_al_2025/Upadated_scripts_after_review_11_2025/Figure1/gibson_barcode_diversity.csv')

# Compute uniformity and compare to Camille CaERG11 library
uniformity_df = compute_uniformity(barcode_mut)
print(uniformity_df)
caERG11 = pd.read_csv('/home/alicia-pageau/Documents/antifungal_project/PDR1/CaERG11_diversity_confirm_CB_2023.csv')
caERG11_uniformity = compute_uniformity(caERG11, count_col='nbr_reads', group_col='Fragment')
print(caERG11_uniformity)
t0_dms = pd.read_csv('/home/alicia-pageau/Documents/gyoza_barcodes_10_2025/results/df/annotated_readcounts/Barc_T0_A_POSA_annot_rc.csv')
t0_dms_uniformity = compute_uniformity(t0_dms, count_col='readcount', group_col='Mutated_seq')
print(t0_dms_uniformity)

barcode_mut["Fragment_num"] = barcode_mut["Fragment"].str.extract(r"F(\d+)").astype(int)
barcode_mut = barcode_mut.sort_values("Fragment_num")
plt.figure(figsize=(10, 4))
ax = sns.boxplot(data=barcode_mut, x='Fragment',y='reads')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_xlabel("Fragment")
ax.set_ylabel("Reads")
plt.show()

#%% Figures (codon)
# Add position mutation in the whole protein seq
barcode_mut['position'] = ((barcode_mut['Fragment'].str.replace("F","").astype(int) -1)*25) + barcode_mut['mut_codon']

# Create a full DataFrame with all combinations
fragments = set(barcode_mut['Fragment'].unique())  # Convert to set for faster lookup
full_index = [(frag, codon) for frag, length in frag_len.items() if frag in fragments for codon in range(1, (length//3) + 1)]
full_df = pd.DataFrame(full_index, columns=['Fragment', 'mut_codon'])
full_df = pd.merge(full_df, df_ref, how = 'left') # this is the same as df_ref...

barcode_mut = full_df.merge(barcode_mut, how='left')
barcode_mut['Fragment'] = pd.Categorical(barcode_mut['Fragment'], categories=natsorted(barcode_mut['Fragment'].unique()), ordered=True)
barcode_mut = barcode_mut.sort_values(['Fragment', 'mut_codon'])


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
#%% Compute barcode diversity (AA)
# Merge reference table with barcode_mut
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
#df.to_csv('/home/alicia-pageau/Documents/antifungal_project/PDR1/gibson/05_analyze/gibson_barcode_diversity.csv')


#%% Figures (AA)
rM_read_depth=pd.read_csv(f"{wkdir}/gibson_all_fragments_sequencing_read_20102025.csv")

barcode_count = rM_read_depth.groupby(['Fragment', 'barcode']).agg(barcode_reads = ('core','nunique')).reset_index()
unique_barcode = barcode_count[barcode_count['barcode_reads'] == 1].drop(columns='barcode_reads')
rM_unique_barcode = pd.merge(rM_read_depth, unique_barcode, on = ['Fragment','barcode'])

# Barcode diversity - Nb unique barcode per AA
rM_unique_barcode['mutation_aa'] = rM_unique_barcode['mutation'].map(codon_to_aa)
rM_unique_barcode['ref_aa'] = rM_unique_barcode['ref_codon'].map(codon_to_aa)

barcode_mut = rM_unique_barcode.groupby(['Fragment','mutation_aa','mut_codon']).agg(barcode_per_mut = ('barcode','nunique')).reset_index()
barcode_mut = barcode_mut[barcode_mut['mutation_aa'] != 'REF']
barcode_mut = barcode_mut.fillna(0)

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
df = df.fillna(0)

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

## PLOT nb barcode per mutation
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
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3.8), constrained_layout=True, dpi=1000)
    axes = axes.flatten() if len(subset) > 1 else [axes]

    for i, fragment in enumerate(subset):
        ax = axes[i]
        ax.set_aspect(1)
        
        # Subset data for the fragment
        df_subset = (barcode_mut[barcode_mut['Fragment'] == fragment]
                     #.assign(barcode_per_mut=lambda d: np.log2(d['barcode_per_mut']))
                     .pivot(index='mutation_aa', columns='position', values='barcode_per_mut'
        ))
        df_subset = df_subset.loc[aa_order]
        
        # Get % of clipped barcodes
        barcode_counts = df_subset.fillna(0).values
        total_barcodes = np.sum(barcode_counts)
        clipped_counts = np.sum(np.maximum(barcode_counts - vmax, 0))
        percent_clipped = (clipped_counts / total_barcodes) * 100
        print(f"Percentage of clipped barcodes for {fragment}: {percent_clipped:.2f}%")
        
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
                ax.add_patch(plt.Rectangle((x_pos + 0.05, y_pos + 0.05), 0.9, 0.9, color='gray', lw=0.1))
    
    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    # Show plots
    plt.show()



