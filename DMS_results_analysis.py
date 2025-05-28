#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 10:57:41 2025

@author: alicia
"""
import pandas as pd
import numpy as np
import seaborn as sns
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import random
import ast

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
barcodes = pd.read_csv("/home/alicia/Documents/gyoza_barcodes/gyoza/results/df/all_scores.csv")
barcodes = barcodes.loc[barcodes['confidence_score'] == 1] # already only replicate A
barcodes = barcodes.groupby(['Drug', 'Fragment', 'Replicate', 'aa_pos', 'alt_aa', 
                 'Nham_aa', 'mutation_type'], as_index=False)['s_T2_T0'].median()

F13 = pd.read_csv("/home/alicia/Documents/gyoza_F13_A/gyoza/results/df/all_scores.csv")
F13 = F13.loc[F13['confidence_score'] == 1] # already only replicate A
F13 = F13.loc[F13['aa_pos'].isin(range(301,326))]

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
fungAMR = pd.read_csv('/home/alicia/Documents/antifungal_project/PDR1/pdr1_F13_res_mutants.csv',sep=';')
fungAMR[['wt_AA', 'position', 'alt_AA']] = fungAMR['mutation'].str.extract(r'(\D+)(\d+)(\D+)')
fungAMR['position'] = fungAMR['position'].astype('Int64')

alignment = pd.read_excel("/home/alicia/Documents/antifungal_project/PDR1/pdr1_alignment_position.xlsx")
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
        df_filtered = df[df['Drug'] == drug].pivot(
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
        
#%% Saturation analysis
# How many barcode are needed per mutation to get good resolution?
barcodes = pd.read_csv("/home/alicia/Documents/antifungal_project/PDR1/gyoza/all_scores_310325.csv")
barcodes = barcodes.loc[barcodes['confidence_score'] == 1] # already only replicate A
# Set barcodes "ID"
barcodes['barcode_ID'] = barcodes.groupby(['Species','Drug','Fragment','Replicate','pos','alt_codons']).cumcount() + 1

F13 = pd.read_csv("/home/alicia/Documents/gyoza_F13_A/gyoza/results/df/all_scores.csv")
F13 = F13.loc[F13['confidence_score'] == 1] # already only replicate A
F13 = F13.loc[F13['aa_pos'].isin(range(301,326))]

barcode = list(range(1,11))
i=10
S = []
for n in range(100):
    random.seed(i)
    sub = random.sample(barcode,10)
    S.append(sub)
    i+=1

for drug in drugs:
    results = []
    barcodes_filtered = barcodes[barcodes['Drug'] == drug]
    
    # Loop through each barcode set in S
    for i, barcode_set in enumerate(S):
        print(f"Processing barcode set {i+1}/{len(S)}")
    
        # Compute correlations for incremental barcode counts
        for nbarcode in range(1, len(barcode_set) + 1):
            current_barcodes = barcode_set[:nbarcode]
            df_subset = barcodes_filtered[barcodes_filtered['barcode_ID'].isin(current_barcodes)]
    
            subset_avg = df_subset.groupby(['Drug', 'Fragment', 'aa_pos', 'alt_aa', 
                                            'Nham_aa', 'mutation_type'], as_index=False)['s_T2_T0'].median()
    
            merged = pd.merge(subset_avg, F13, 
                              on=['Drug', 'Fragment', 'aa_pos', 'alt_aa', 'Nham_aa', 'mutation_type'])
    
            if len(merged) > 1:
                corr, _ = pearsonr(merged['s_T2_T0_x'], merged['s_T2_T0_y'])
            else:
                corr = np.nan
    
            results.append({
                'barcode_set_id': i,
                'n_barcodes': nbarcode,
                'correlation': corr
            })
    
    # Convert results to DataFrame
    corr_df = pd.DataFrame(results)    
    
    plt.figure(figsize=(6, 4.5), dpi = 300)
    sns.pointplot(
        data=corr_df,
        x='n_barcodes',
        y='correlation',
        errorbar='ci',
        color='black',
        marker='o',
        join=True
    )
    plt.ylim(0, 1)
    plt.xlabel('Number of barcodes', fontsize = 12)
    plt.ylabel('Correlation coefficient\n(DMS F13 sequence Vs DMS DNA barcode)', fontsize = 12)
    plt.title(f'{drug}', fontsize = 12)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


