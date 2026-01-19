#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 10:57:41 2025

@author: alicia-pageau
10-11-2025
PDR1 F13 DMS Analysis
"""
import pandas as pd
import numpy as np
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from itertools import combinations
import random
import ast
import os

wkdir = '/home/alicia-pageau/Documents/antifungal_project/PDR1'
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


#%% Heatmaps from gyoza all_score.csv result
barcodes = pd.read_csv("/home/alicia-pageau/Documents/gyoza_barcodes_10_2025/results/df/all_scores.csv")
barcodes = barcodes.loc[barcodes['confidence_score'] == 1] 
barcodes = barcodes.loc[barcodes['aa_pos'] != 'not-applicable']
barcodes['aa_pos'] = barcodes['aa_pos'].astype(int)
barcodes = barcodes.loc[barcodes['aa_pos'].isin(range(301,326))]

barcodes = barcodes.groupby(['Drug', 'aa_pos', 'alt_aa', 
                 'Nham_aa', 'mutation_type'], as_index=False)['s_T2_T0'].median()

F13 = pd.read_csv("/home/alicia-pageau/Documents/gyoza_F13_10_2025/results/df/all_scores.csv")
F13 = F13.loc[F13['confidence_score'] == 1] 
F13 = F13.loc[F13['aa_pos'] != 'not-applicable']
F13['aa_pos'] = F13['aa_pos'].astype(int)
F13 = F13.loc[F13['aa_pos'].isin(range(301,326))]
F13 = F13[F13['Replicate']=='A']

F13 = F13.groupby(['Drug', 'aa_pos', 'alt_aa', 
                 'Nham_aa', 'mutation_type'], as_index=False)['s_T2_T0'].median()

# Create reference DataFrame
dref = list(SeqIO.parse("references_PDR1_corrected.fasta","fasta"))

R = {rec.id.split('_')[1]: str(rec.seq).upper() for rec in dref}
aa_data = []
for seq_id, sequence in R.items():
    for i in range(0, len(sequence) - 2, 3):  # Iterate in steps of 3
        codon = sequence[i:i+3]  # Extract codon
        codon_position = i // 3 + 1  # Convert base position to codon number
        
        # 🔹 Convert codon to amino acid 
        amino_acid = codon_to_aa.get(codon)
        
        aa_data.append([seq_id, codon_position, codon, amino_acid])

df_ref = pd.DataFrame(aa_data, columns=["Fragment", "mut_codon", "ref_codon", "ref_aa"])
df_ref['position'] = ((df_ref['Fragment'].str.replace("F","").astype(int) -1)*25) + df_ref['mut_codon']

# FungAMR mutations
fungAMR = pd.read_csv('/home/alicia-pageau/Documents/antifungal_project/PDR1/pdr1_F13_res_mutants.csv',sep=';')
fungAMR[['wt_AA', 'position', 'alt_AA']] = fungAMR['mutation'].str.extract(r'(\D+)(\d+)(\D+)')
fungAMR['position'] = fungAMR['position'].astype('Int64')

alignment = pd.read_excel("/home/alicia-pageau/Documents/antifungal_project/PDR1/pdr1_alignment_position.xlsx")
alignment = alignment.rename(columns={
    'Saccharomyces cerevisiae_position': 'scer_pos',
    'Saccharomyces cerevisiae_wt_AA': 'scer_aa',
    'Saccharomyces paradoxus_position': 'spar_pos',
    'Saccharomyces paradoxus_wt_AA': 'spar_aa',
    'Nakaseomyces glabratus_position': 'ngla_pos',
    'Nakaseomyces glabratus_wt_AA': 'ngla_aa'
})
alignment['spar_pos'] = alignment['spar_pos'].astype('Int64')
alignment['ngla_pos'] = alignment['ngla_pos'].astype('Int64')
alignment['scer_pos'] = alignment['scer_pos'].astype('Int64')

def map_to_cerevisiae(mutations_df, alignment_df):
    # Prepare output list
    mapped_data = []

    for _, row in mutations_df.iterrows():
        species = row['species']
        mut_pos = row['position']

        if 'paradoxus' in species:
            match = alignment[alignment['spar_pos'] == mut_pos]
        elif 'glabrata' in species:
            match = alignment[alignment['ngla_pos'] == mut_pos]
        elif 'cerevisiae' in species:
            match = alignment[alignment['scer_pos'] == mut_pos]
        else:
            match = pd.DataFrame()  # Unknown species

        if not match.empty:
            scer_pos = match['scer_pos'].values[0]
            scer_aa = match['scer_aa'].values[0]
        else:
            scer_pos = None
            scer_aa = None

        mapped_data.append({
            **row,
            'scer_position': scer_pos,
            'scer_wt_aa': scer_aa
        })

    return pd.DataFrame(mapped_data)

fungAMR['species'] = fungAMR['species'].apply(ast.literal_eval)
fungAMR['species'] = fungAMR['species'].str[0]  # if species is in a list format
fungAMR = map_to_cerevisiae(fungAMR, alignment)
fungAMR['scer_position'] = fungAMR['scer_position'].astype('Int64')
fungAMR = fungAMR.loc[(fungAMR['scer_position']>=301)&(fungAMR['scer_position']<=325)]

# PLOT
drugs = ['CTL', 'NQO', 'POSA']
dfs = [("barcodes", barcodes), ("F13", F13)]

for df_name, df in dfs:
    for drug in drugs:
        # Filter the DataFrame for the current drug
        df_filtered = df[(df['Drug'] == drug)].pivot(
            index='alt_aa', columns='aa_pos', values='s_T2_T0'
        )
        df_filtered = df_filtered.loc[aa_order]
    
        # Plot heatmap s_T2_T0
        fig, ax = plt.subplots(figsize=(8, 4.8), dpi = 300)
        sns.heatmap(df_filtered, cmap="coolwarm",linewidths=0.5,
                    vmin = -0.6,vmax=0.6)
        
        # ✅ Add gray squares for WT codons
        for idx, row in df_ref[df_ref['Fragment'] == 'F13'].iterrows():
            ref_codon = str(row['ref_aa'])  # Ensure it's a string
            mut_codon = row['position']  # This is a position (number)

            y_pos = df_filtered.index.get_loc(ref_codon)  # Get row index
            x_pos = df_filtered.columns.get_loc(mut_codon)  # Get column index
                
            # ✅ Add gray square for WT codons
            #ax.add_patch(plt.Rectangle((x_pos + 0.05, y_pos +0.05), 0.95, 0.95, color='black', lw=0.1))
            ax.plot(x_pos + 0.5, y_pos +0.5, 'o', color='black', markersize=5)
                
        
        if drug == 'POSA':
            for _, row in fungAMR.iterrows():
                species = row['species']
                
                alt_codon = str(row['alt_AA'])  # Ensure it's a string
                mut_codon = row['scer_position']  # This is a position (number)
    
                y_pos = df_filtered.index.get_loc(alt_codon)  # Get row index
                x_pos = df_filtered.columns.get_loc(mut_codon)  # Get column index
                
                if species == 'Nakaseomyces glabrata':
                    ax.add_patch(plt.Rectangle((x_pos, y_pos), 1, 1, facecolor='none', edgecolor='purple', lw=2.5))
                elif mut_codon == 306 and alt_codon in ['C', 'H']:
                    ax.add_patch(plt.Rectangle((x_pos+0.1, y_pos+0.1), 0.8, 0.8, facecolor='none', edgecolor='gray', lw=3))
                else:
                    ax.add_patch(plt.Rectangle((x_pos, y_pos), 1, 1, facecolor='none', edgecolor='gray', lw=2.5))
        
        cbar = ax.collections[0].colorbar
        cbar.ax.set_ylabel('Selection coefficient', rotation=270, labelpad=15, fontsize = 12)
        ax.set_title(f"DMS {df_name} for {drug} (s_T2_T0)")
        ax.set_xlabel("Position", fontsize = 12)
        ax.set_ylabel("")
        plt.show()
        
#%% Correlation F13 VS barcode all usable barcode
barcodes = pd.read_csv("/home/alicia-pageau/Documents/gyoza_barcodes_10_2025/results/df/all_scores.csv")
barcodes = barcodes.loc[barcodes['confidence_score'] == 1] # already only replicate A
barcodes = barcodes.loc[barcodes['aa_pos'] != 'not-applicable']
barcodes['aa_pos'] = barcodes['aa_pos'].astype(int)
barcodes = barcodes.loc[barcodes['aa_pos'].isin(range(301,326))]
barcodes['nbarcodes'] = barcodes.groupby(['Drug','aa_pos','alt_aa'])['barcode'].transform('nunique')
barcodes.rename(columns={'s_T2_T0':'s_T2_T0_barcodes'},inplace = True)

F13 = pd.read_csv("/home/alicia-pageau/Documents/gyoza_F13_10_2025/results/df/all_scores.csv")
F13 = F13.loc[F13['confidence_score'] == 1] 
F13 = F13.loc[F13['aa_pos'] != 'not-applicable']
F13['aa_pos'] = F13['aa_pos'].astype(int)
F13 = F13.loc[F13['aa_pos'].isin(range(301,326))]
F13.rename(columns={'s_T2_T0':'s_T2_T0_F13'},inplace = True)

# Correlation per AA
barcodes_aa = barcodes.groupby(['Drug', 'aa_pos', 'alt_aa', 
                 'Nham_aa', 'mutation_type', 'nbarcodes'], as_index=False)['s_T2_T0_barcodes'].median()
F13_aa = F13.groupby(['Drug', 'aa_pos', 'alt_aa', 
                 'Nham_aa', 'mutation_type'], as_index=False)['s_T2_T0_F13'].median()

data = pd.merge(barcodes_aa,F13_aa, on = ['Drug', 'aa_pos', 'alt_aa', 'Nham_aa', 'mutation_type'])
data['s_diff'] = abs(data['s_T2_T0_barcodes']-data['s_T2_T0_F13'])
data = data[data['Drug']=='POSA']
corr = []

# Make a scatter plot between experiments
plt.figure(figsize=(6,6), dpi=300)
ax = sns.scatterplot(
    data=data,
    x='s_T2_T0_barcodes',   # first experiment
    y='s_T2_T0_F13',   # second experiment
    hue="nbarcodes",
    s=60,
    alpha=0.8,
)
for i, drug in enumerate(data["Drug"].unique()):
    subset = data[data["Drug"] == drug]
    # make sure you pass Series (1D), not DataFrame (2D)
    x = subset["s_T2_T0_barcodes"].dropna()
    y = subset["s_T2_T0_F13"].dropna()

    if len(x) > 1 and len(y) > 1:  # need at least 2 points
        r, p = pearsonr(x, y)
        corr.append({
            "Drug": drug,
            "pearson_r": r,
            "p_value": p,
            "n_points": len(subset)
        })
        plt.text(
            0.05,
            0.95 - i*0.05,
            f"r = {r:.2f}",
            transform=plt.gca().transAxes,
            fontsize=10
        )

plt.xlabel("Selection coefficient (DNA barcodes)")
plt.ylabel("Selection coefficient (PDR1 F13 sequence)")
plt.axis("equal")
plt.legend(title="Number of barcode\nper mutation")
plt.tight_layout()
plt.show()

corr_df = pd.DataFrame(corr)

#%% Internal reproducibility (codon and barcodes)
barcodes = pd.read_csv("/home/alicia-pageau/Documents/gyoza_barcodes_10_2025/results/df/all_scores.csv")
barcodes = barcodes.loc[barcodes['confidence_score'] == 1] 
barcodes = barcodes.loc[barcodes['aa_pos'] != 'not-applicable']
barcodes['aa_pos'] = barcodes['aa_pos'].astype(int)
barcodes = barcodes[barcodes['Drug']=='POSA']
barcodes['barcode_ID'] = barcodes.groupby(['aa_pos','mutation_alt_codons']).cumcount() + 1

barcodes['pattern'] = barcodes['barcode'].apply(lambda x: x[2:14] + "_" + x[16:28])
barcodes['pattern_ID'] = barcodes.groupby(['pattern']).ngroup() + 1


F13 = pd.read_csv("/home/alicia-pageau/Documents/gyoza_F13_10_2025/results/df/all_scores.csv")
F13 = F13.loc[F13['confidence_score'] == 1] 
F13 = F13.loc[F13['aa_pos'] != 'not-applicable']
F13['aa_pos'] = F13['aa_pos'].astype(int)
F13 = F13.loc[F13['aa_pos'].isin(range(301,326))]
F13 = F13[F13['Drug']=='POSA']

def corrfunc(x, y, **kws):
    mask = ~np.isnan(x) & ~np.isnan(y)
    r, _ = pearsonr(x[mask], y[mask])
    ax = plt.gca()
    ax.annotate(f'{r:.2f}', xy=(.5, .5), xycoords='axes fraction', fontsize=16)
    ax.set_xlim(-0.5, 0.65)
    ax.set_ylim(-0.5, 0.65)
    ax.set_axis_off()
    
def scatterfunc(x, y, **kws):
    if x.name == y.name:
        return
    ax = plt.gca()
    sns.scatterplot(x=x, y=y, **kws)
    ax.set_xlim(-0.5, 0.65)
    ax.set_ylim(-0.5, 0.65)
    
def label_diag(x, **kwargs):
    ax = plt.gca()
    idx = np.where(np.array(g.diag_axes) == ax)[0][0]
    label = g.x_vars[idx]
    ax.text(0.5, 0.5, label, ha='center', va='center', transform=ax.transAxes,
            fontsize=18, fontweight='bold')
    ax.set_axis_off()

def pairwise_corr_with_n(df, columns=None):
    """
    Compute pairwise Pearson correlation and sample size (n) between specified columns.
    
    Parameters:
        df : pd.DataFrame
            Input DataFrame with numeric columns.
        columns : list or None
            List of columns to compare. If None, all numeric columns are used.
    
    Returns:
        pd.DataFrame with columns ["col1", "col2", "n", "r"]
    """
    if columns is None:
        columns = df.select_dtypes(include='number').columns.tolist()
        
    results = []
    
    for i, c1 in enumerate(columns):
        for j, c2 in enumerate(columns):
            if j <= i:
                continue  # avoid duplicate pairs and self-comparison
            x = df[c1]
            y = df[c2]
            valid = x.notna() & y.notna()
            n = valid.sum()
            r = x[valid].corr(y[valid]) if n > 1 else None
            results.append((c1, c2, n, r))
    
    return pd.DataFrame(results, columns=["col1", "col2", "n", "r"])


########## Barcodes ##########
# Within-codon reproducibility (barcodes)
codon_stats = (
    barcodes.groupby(['aa_pos','mutation_alt_codons'])['s_T2_T0']
      .agg(['mean','median','std', 'count','min','max'])
      .assign(CV=lambda x: x['std'] / x['mean'])
).reset_index()

pivot = barcodes.pivot_table(index=['aa_pos','mutation_alt_codons'], columns='barcode_ID', values='s_T2_T0')
pivot = pivot.loc[:, pivot.nunique(dropna=True) > 2]
pivot = pivot.add_prefix("barcode ")

g = sns.PairGrid(pivot.iloc[:, :5])
g.map_lower(scatterfunc)
g.map_diag(label_diag)
g.map_upper(corrfunc)
g.set(xlabel=None, ylabel=None)
plt.suptitle("Barcode reproducibility (DNA barcode)", y=1.02)
plt.show()

pair_corr = pairwise_corr_with_n(pivot)
mean_r = pair_corr['r'].mean()
median_r = pair_corr['r'].median()
q1 = pair_corr['r'].quantile(0.25)
q3 = pair_corr['r'].quantile(0.75)
iqr = q3 - q1
print(f"Mean: {mean_r:.3f}")
print(f"Median: {median_r:.3f}")
print(f"IQR: {iqr:.3f} (Q1={q1:.3f}, Q3={q3:.3f})")

# Pattern reproducibility
pattern_median = (
    barcodes.groupby(['aa_pos','mutation_alt_aa','mutation_alt_codons','pattern'])['s_T2_T0']
      .median().reset_index()
)
pattern_median['pattern_ID'] = pattern_median.groupby(['aa_pos','mutation_alt_aa','mutation_alt_codons']).cumcount() + 1

pivot = pattern_median.pivot_table(index=['aa_pos','mutation_alt_aa','mutation_alt_codons'], columns='pattern_ID', values='s_T2_T0')
pivot = pivot.loc[:, pivot.nunique(dropna=True) > 2]
pivot = pivot.add_prefix("pattern ")

pair_corr = pairwise_corr_with_n(pivot)
mean_r = pair_corr['r'].mean()
median_r = pair_corr['r'].median()
q1 = pair_corr['r'].quantile(0.25)
q3 = pair_corr['r'].quantile(0.75)
iqr = q3 - q1
print(f"Mean: {mean_r:.3f}")
print(f"Median: {median_r:.3f}")
print(f"IQR: {iqr:.3f} (Q1={q1:.3f}, Q3={q3:.3f})")

g = sns.PairGrid(pivot.iloc[:,:5])
g.map_lower(scatterfunc)
g.map_diag(label_diag)
g.map_upper(corrfunc)
g.set(xlabel=None, ylabel=None)
plt.suptitle("Pattern reproducibility (DNA barcode)", y=1.02)
plt.show()


# Within-aa reproducibility (Synonymous codon)
aa_median = (
    barcodes.groupby(['aa_pos','mutation_alt_aa','mutation_alt_codons'])['s_T2_T0']
      .median().reset_index()
)
aa_median['codon_ID'] = aa_median.groupby(['aa_pos','mutation_alt_aa']).cumcount() + 1

pivot = aa_median.pivot_table(index=['aa_pos','mutation_alt_aa'], columns='codon_ID', values='s_T2_T0')
pivot = pivot.loc[:, pivot.nunique(dropna=True) > 2]
pivot = pivot.add_prefix("codon ")

pair_corr = pairwise_corr_with_n(pivot)
mean_r = pair_corr['r'].mean()
median_r = pair_corr['r'].median()
q1 = pair_corr['r'].quantile(0.25)
q3 = pair_corr['r'].quantile(0.75)
iqr = q3 - q1
print(f"Mean: {mean_r:.3f}")
print(f"Median: {median_r:.3f}")
print(f"IQR: {iqr:.3f} (Q1={q1:.3f}, Q3={q3:.3f})")

g = sns.PairGrid(pivot)
g.map_lower(scatterfunc)
g.map_diag(label_diag)
g.map_upper(corrfunc)
g.set(xlabel=None, ylabel=None)
plt.suptitle("Synonymous codon reproducibility (DNA barcode)", y=1.02)
plt.show()


######### F13 ##########
# Within-aa reproducibility (Synonymous codon)
aa_median = (
    F13.groupby(['aa_pos','mutation_alt_aa','mutation_alt_codons'])['s_T2_T0']
      .median().reset_index()
)
aa_median['codon_ID'] = aa_median.groupby(['aa_pos','mutation_alt_aa']).cumcount() + 1

pivot = aa_median.pivot_table(index=['aa_pos','mutation_alt_aa'], columns='codon_ID', values='s_T2_T0')
pivot = pivot.loc[:, pivot.nunique(dropna=True) > 2]
pivot = pivot.add_prefix("codon ")

pair_corr = pairwise_corr_with_n(pivot)
mean_r = pair_corr['r'].mean()
median_r = pair_corr['r'].median()
q1 = pair_corr['r'].quantile(0.25)
q3 = pair_corr['r'].quantile(0.75)
iqr = q3 - q1
print(f"Mean: {mean_r:.3f}")
print(f"Median: {median_r:.3f}")
print(f"IQR: {iqr:.3f} (Q1={q1:.3f}, Q3={q3:.3f})")

g = sns.PairGrid(pivot)
g.map_lower(scatterfunc)
g.map_diag(label_diag)
g.map_upper(corrfunc)
g.set(xlabel=None, ylabel=None)
plt.suptitle("Synonymous codon reproducibility (F13 sequence)", y=1.02)
plt.show()

# Replicate reproducibility
rep_median = (
    F13.groupby(['aa_pos', 'mutation_alt_aa', 'mutation_alt_codons', 'Replicate'])['s_T2_T0']
       .median()
       .reset_index()
)

pivot = rep_median.pivot_table(index=['aa_pos','mutation_alt_aa'], columns='Replicate', values='s_T2_T0')
g = sns.PairGrid(pivot, diag_sharey=True)
g.map_lower(scatterfunc)
g.map_diag(label_diag)
g.map_upper(corrfunc)
g.set(xlabel=None, ylabel=None)
plt.suptitle("Replicates reproducibility (F13 sequence)", y=1.02)
plt.show()

pair_corr = pairwise_corr_with_n(pivot)
mean_r = pair_corr['r'].mean()
median_r = pair_corr['r'].median()
q1 = pair_corr['r'].quantile(0.25)
q3 = pair_corr['r'].quantile(0.75)
iqr = q3 - q1
print(f"Mean: {mean_r:.3f}")
print(f"Median: {median_r:.3f}")
print(f"IQR: {iqr:.3f} (Q1={q1:.3f}, Q3={q3:.3f})")


#%% Selection coeff STOP codon 
barcodes = pd.read_csv("/home/alicia-pageau/Documents/gyoza_barcodes_10_2025/results/df/all_scores.csv")
barcodes = barcodes.loc[barcodes['confidence_score'] == 1] # already only replicate A
barcodes = barcodes.loc[barcodes['aa_pos'] != 'not-applicable']
barcodes['aa_pos'] = barcodes['aa_pos'].astype(int)
barcodes = barcodes.loc[barcodes['aa_pos'].isin(range(301,326))]
barcodes = barcodes.groupby(['Drug', 'aa_pos', 'alt_aa', 
                 'Nham_aa', 'mutation_type'], as_index=False)['s_T2_T0'].median()
barcodes_stop_codon = barcodes[barcodes['alt_aa']=='*']
barcodes_stop_codon['experiment'] = 'barcode'
barcodes_stop_codon['Replicate'] = 'A'


F13 = pd.read_csv("/home/alicia-pageau/Documents/gyoza_F13_10_2025/results/df/all_scores.csv")
F13 = F13.loc[F13['confidence_score'] == 1] 
F13 = F13.loc[F13['aa_pos'] != 'not-applicable']
F13['aa_pos'] = F13['aa_pos'].astype(int)
F13 = F13.loc[F13['aa_pos'].isin(range(301,326))]
F13_stop_codon = F13[F13['alt_aa']=='*']
F13_stop_codon['experiment'] = 'F13'

stop_codon = pd.concat([barcodes_stop_codon,F13_stop_codon])
stop_codon = stop_codon[stop_codon['Drug']=='POSA']


plt.figure(figsize=(4,5), dpi = 300)

sns.stripplot(
    data=stop_codon,
    x="experiment",
    y="s_T2_T0",
    hue="Replicate",
    dodge=False,     # separate points per hue
    jitter=True,    # spread points for visibility
    alpha=1,      # transparency
    size=4,         # point size
    linewidth=0.5,
    edgecolor="gray"
)

sns.violinplot(
    data=stop_codon,
    x="experiment",
    y="s_T2_T0",
    inner=None,     # remove inner bars since we add points
    #cut=0,          # avoid extending beyond data range
    linewidth=1,
    alpha=0.5,
    color="gray"
)

handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:len(stop_codon["Replicate"].unique())],
           labels[:len(stop_codon["Replicate"].unique())],
           title="Replicate")
plt.xticks(ticks=[0, 1], labels=["DNA barcodes", "PDR1 F13 sequence"])
plt.xlabel("")
plt.ylabel("Selection coefficient")
plt.tight_layout()
plt.show()

    

#%% Saturation analysis (How many barcode are needed per mutation to get good resolution?)
barcodes = pd.read_csv("/home/alicia-pageau/Documents/gyoza_barcodes_10_2025/results/df/all_scores.csv")
barcodes = barcodes.loc[barcodes['confidence_score'] == 1] # already only replicate A
barcodes = barcodes.loc[barcodes['aa_pos'] != 'not-applicable']
barcodes['aa_pos'] = barcodes['aa_pos'].astype(int)
barcodes = barcodes.loc[barcodes['aa_pos'].isin(range(301,326))]
barcodes['barcode_ID'] = barcodes.groupby(['Drug','aa_pos','alt_aa']).cumcount() + 1
barcodes['nbarcodes'] = barcodes.groupby(['Drug','aa_pos','alt_aa'])['barcode'].transform('nunique')

F13 = pd.read_csv("/home/alicia-pageau/Documents/gyoza_F13_10_2025/results/df/all_scores.csv")
F13 = F13.loc[F13['confidence_score'] == 1] 
F13 = F13.loc[F13['aa_pos'] != 'not-applicable']
F13['aa_pos'] = F13['aa_pos'].astype(int)
F13 = F13.loc[F13['aa_pos'].isin(range(301,326))]
F13 = F13.groupby(['Drug', 'aa_pos', 'alt_aa','Nham_aa', 'mutation_type'], 
                  as_index=False)['s_T2_T0'].median()


def sample_with_reuse(values, k):
    """
    Always include all available values, then sample extra (with replacement)
    until k is reached.
    """
    values = list(values)
    n = len(values)
    if n >= k:
        return np.random.choice(values, size=k, replace=False)
    else:
        sampled = values.copy()
        extra = np.random.choice(values, size=k-n, replace=True)
        sampled.extend(extra)
        return sampled

def bootstrap_correlation(df_barcode, df_dms, k, barcode_tot=None, drug='POSA', n_iter=100, method="pearson"):
    """
    Subsample k barcodes per mutant using the 'always include all' rule.
    """
    results = []
    df_barcode = df_barcode[df_barcode['Drug'] == drug]
    if barcode_tot is not None:
        df_barcode = df_barcode[df_barcode['nbarcodes'] == barcode_tot]
    df_dms = df_dms[df_dms['Drug'] == drug]
    
    for i in range(n_iter):
        pooled = []
        #picked = []
        for (aa_pos, alt_aa), sub in df_barcode.groupby(['aa_pos', 'alt_aa']):
            barcodes = sub["barcode_ID"].values
            sampled_barcodes = sample_with_reuse(barcodes, k)
            sampled_coeffs = sub.set_index("barcode_ID").loc[sampled_barcodes, "s_T2_T0"].values
            pooled.append([aa_pos, alt_aa, np.median(sampled_coeffs)])
        
        pooled_df = pd.DataFrame(pooled, columns=['aa_pos', 'alt_aa', "s_median"])
        
        merged = pooled_df.merge(df_dms, on=['aa_pos', 'alt_aa'])
        n_variants = merged[['aa_pos', 'alt_aa']].drop_duplicates().shape[0]

        corr = merged["s_median"].corr(merged["s_T2_T0"], method=method)

        results.append({
            'iter': i,
            'barcode_tot' : barcode_tot,
            'n_barcodes': k,
            'correlation': corr,
            'n_variants' : n_variants
        })
    
    return np.array(results)

all_results = []
for k in range(1, 11):
    res = bootstrap_correlation(barcodes, F13, k, n_iter=100)
    all_results.extend(res)  # append all dicts
corr_df = pd.DataFrame(all_results)

plt.figure(figsize=(9, 6.75), dpi = 300)
sns.boxplot(
    data=corr_df,
    x='n_barcodes',
    y='correlation',
    color='lightgray'
)
plt.ylim(0, 1)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel('Number of barcodes', fontsize = 20)
plt.ylabel('Correlation coefficient\n(DMS F13 sequence Vs DMS barcode)', fontsize = 20)
plt.grid(True)
plt.tight_layout()
plt.show()

# Check if barcode diversity affects quality
barcode_tot_list = [5, 10, 15]
all_results = []
for barcode_tot in barcode_tot_list:
    for k in range(1, 16):
        res = bootstrap_correlation(barcodes, F13, k, 
                                    barcode_tot=barcode_tot, n_iter=100)
        all_results.extend(res)
corr_df_all = pd.DataFrame(all_results)
corr_df_all.loc[corr_df_all['n_barcodes'] > corr_df_all['barcode_tot'],'correlation'] = np.nan

plt.figure(figsize=(9, 6.75), dpi=300)
ax =sns.pointplot(
    data=corr_df_all,
    x='n_barcodes',
    y='correlation',
    hue='barcode_tot',  # different lines for different barcode_tot
    errorbar='sd',
    join=True
)
plt.ylim(0, 1)
plt.xlabel('Number of barcodes', fontsize=14)
plt.ylabel('Correlation coefficient\n(DMS F13 sequence Vs DMS DNA barcode)', fontsize=14)
plt.grid(True)
hue_colors = {}
for lh in ax.get_legend_handles_labels()[0]:  # handles in legend
    hue = lh.get_label()
    hue_colors[hue] = lh.get_color()
for i, barcode_tot in enumerate(barcode_tot_list):
    n_variants = corr_df_all[corr_df_all['barcode_tot'] == barcode_tot]['n_variants'].iloc[0]
    ax.text(
        x=0 + i*2,
        y=1.02,
        s=f'{barcode_tot}: n={n_variants}',
        color=hue_colors[str(barcode_tot)],  # convert to str if needed
        fontsize=10
    )
plt.tight_layout()
plt.show()

