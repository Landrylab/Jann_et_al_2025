#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 16:09:52 2026

@author: alicia-pageau
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import os

wkdir = '/home/alicia-pageau/Documents/antifungal_project/PDR1/00_scripts/Jann_et_al_2025/Upadated_scripts_after_review_11_2025/Figure2/'
os.chdir(wkdir)

# aa order per properties for heatmaps
aa_order = ["*", "P", "G", "C", "Q", "N", "T", "S", "E", "D",
            "K", "H", "R", "W", "Y", "F", "M", "L", "I", "V", "A"]

# %% Import files
df_ref = pd.read_csv('ref.csv', index_col=0)

# DMS gyoza results
barcodes = pd.read_csv("gyoza_barcodes_select_coeff_median.csv", index_col=0)
F13 = pd.read_csv("gyoza_F13_select_coeff_median.csv", index_col=0)

# FungAMR mutations
fungAMR = pd.read_csv('fungAMR_pdr1_F13_res_mutants.csv', index_col=0)

# Barcode saturation
barcode_sat = pd.read_csv('barcode_saturation_analysis.csv', index_col=0)

# %% Figure 2B - Heatmaps from DMS results
drugs = ['POSA']
dfs = [("barcodes", barcodes), ("F13", F13)]

for df_name, df in dfs:
    for drug in drugs:
        df_filtered = df[(df['Drug'] == drug)].pivot(
            index='alt_aa', columns='aa_pos', values='s_T2_T0'
        )
        df_filtered = df_filtered.loc[aa_order]
    
        # Plot heatmap selection coefficient
        fig, ax = plt.subplots(figsize=(8, 4.8), dpi = 300)
        sns.heatmap(df_filtered, cmap="coolwarm",linewidths=0.5,
                    vmin = -0.6,vmax=0.6)
        
        # Add black dot for WT aa
        for idx, row in df_ref[df_ref['Fragment'] == 'F13'].iterrows():
            ref_aa = str(row['ref_aa'])
            mut_aa = row['position']
            
            y_pos = df_filtered.index.get_loc(ref_aa)
            x_pos = df_filtered.columns.get_loc(mut_aa)
            
            ax.plot(x_pos + 0.5, y_pos +0.5, 'o', color='black', markersize=5)
        
        # Highlight orthologous mutation from fungAMR
        if drug == 'POSA':
            for _, row in fungAMR.iterrows():
                species = row['species']
                
                alt_aa = str(row['alt_AA'])
                mut_aa = row['scer_position'] 
    
                y_pos = df_filtered.index.get_loc(alt_aa)
                x_pos = df_filtered.columns.get_loc(mut_aa)
                
                if species == 'Nakaseomyces glabrata':
                    ax.add_patch(plt.Rectangle((x_pos, y_pos), 1, 1, facecolor='none', edgecolor='purple', lw=2.5))
                elif mut_aa == 306 and alt_aa in ['C', 'H']:
                    ax.add_patch(plt.Rectangle((x_pos+0.1, y_pos+0.1), 0.8, 0.8, facecolor='none', edgecolor='gray', lw=3))
                else:
                    ax.add_patch(plt.Rectangle((x_pos, y_pos), 1, 1, facecolor='none', edgecolor='gray', lw=2.5))
        
        cbar = ax.collections[0].colorbar
        cbar.ax.set_ylabel('Selection coefficient', rotation=270, labelpad=15, fontsize = 12)
        ax.set_title(f"DMS {df_name} for {drug}")
        ax.set_xlabel("Position", fontsize = 12)
        ax.set_ylabel("")
        
        # Save figures
        fig.savefig(f"figure_2B_{df_name}",bbox_inches="tight", dpi=300)
        plt.close(fig)

# %% Figure 2C - Barcodes saturation analysis
plt.figure(figsize=(9, 6.75), dpi = 300)
sns.boxplot(
    data=barcode_sat,
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
plt.savefig("figure_2C.png", dpi=300, bbox_inches="tight")
