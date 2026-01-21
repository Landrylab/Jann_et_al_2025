#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:42:54 2026

@author: alicia-pageau
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import os

wkdir = '/home/alicia-pageau/Documents/antifungal_project/PDR1/00_scripts/Jann_et_al_2025/Upadated_scripts_after_review_11_2025/Supplementary_figures/'
os.chdir(wkdir)

# aa order per properties for heatmaps
aa_order = ["*", "P", "G", "C", "Q", "N", "T", "S", "E", "D",
            "K", "H", "R", "W", "Y", "F", "M", "L", "I", "V", "A"]


# %% Supp figure 13 - Correlation all-usable barcodes
F13 = pd.read_csv("S5_data.csv")
F13 = F13.loc[F13['confidence_score'] == 1] 
F13 = F13.loc[F13['aa_pos'] != 'not-applicable']
F13['aa_pos'] = F13['aa_pos'].astype(int)
F13 = F13.loc[F13['aa_pos'].isin(range(301,326))]
F13 = F13[F13['Replicate']=='A']
F13 = F13.groupby(['Drug', 'aa_pos', 'alt_aa', 'Nham_aa', 'mutation_type'], as_index=False)['s_T2_T0'].median()

barcodes = pd.read_csv("S6_data.csv",)
barcodes = barcodes.loc[barcodes['confidence_score'] == 1] 
barcodes = barcodes.loc[barcodes['aa_pos'] != 'not-applicable']
barcodes['aa_pos'] = barcodes['aa_pos'].astype(int)
barcodes = barcodes.loc[barcodes['aa_pos'].isin(range(301,326))]
barcodes['nbarcodes'] = barcodes.groupby(['Drug','aa_pos','alt_aa'])['barcode'].transform('nunique')
barcodes = barcodes.groupby(['Drug', 'aa_pos', 'alt_aa', 'Nham_aa', 'mutation_type', 'nbarcodes'], as_index=False)['s_T2_T0'].median()

data = pd.merge(barcodes,F13, on = ['Drug', 'aa_pos', 'alt_aa', 'Nham_aa', 'mutation_type'])
data = data[data['Drug']=='POSA']

# Compute correlation
x = data['s_T2_T0_x']
y = data['s_T2_T0_y']
r, p = pearsonr(x, y)

# Make a scatter plot 
plt.figure(figsize=(6,6), dpi=300)
ax = sns.scatterplot(
    data=data,
    x='s_T2_T0_x',
    y='s_T2_T0_y',
    hue="nbarcodes",
    s=60,
    alpha=0.8,
)
plt.text(
    0.05,
    0.95,
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
#plt.savefig("supp_figure_13.png", dpi=300, bbox_inches="tight")

# %% Supp figure 14 - Barcodes saturation analysis diversity impact
barcode_div_sat = pd.read_csv('S9_data.csv')

plt.figure(figsize=(9, 6.75), dpi=300)
ax = sns.pointplot(
    data=barcode_div_sat,
    x='n_barcodes',
    y='correlation',
    hue='barcode_tot',
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
for i, barcode_tot in enumerate(barcode_div_sat['barcode_tot'].unique().tolist()):
    n_variants = barcode_div_sat[barcode_div_sat['barcode_tot'] == barcode_tot]['n_variants'].iloc[0]
    ax.text(
        x=0 + i*2,
        y=1.02,
        s=f'{barcode_tot}: n={n_variants}',
        color=hue_colors[str(barcode_tot)],
        fontsize=10
    )
plt.tight_layout()
plt.savefig("supp_figure_14.png", dpi=300, bbox_inches="tight")
