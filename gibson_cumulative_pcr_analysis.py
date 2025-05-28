#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from Bio import SeqIO
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
import random

wkdir = '/home/alicia/Documents/antifungal_project/PDR1'
os.chdir(wkdir)

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

def calc_nunique_refs(blist):
    refs = [ele for ele in blist if 'REF' in ele]
    return(len(set(refs)))

fnames = glob.glob(f"{wkdir}/gibson/05_analyze/F*_*.csv")
barcode_cov = 1

# Main table gibson reads
M = []
for fname in fnames:
    df_ = pd.read_csv(fname, sep=",", header=0)
    dbarcore = df_.groupby(['Fragment','pcr','barcode','core']).agg(barcore_reads = ('seq','count')).reset_index()
    df_2 = pd.merge(df_, dbarcore, on = ['Fragment','pcr','barcode','core'], how = 'left')
    dsub = df_2[df_2['barcore_reads'] > barcode_cov].reset_index(drop=True)

    M.append(dsub)

dM = pd.concat(M, ignore_index=True)
dM.to_csv("02_combine_results_"+str(barcode_cov)+".csv", sep=",", header=True, index=False)

df = dM.copy()
df['mutation_aa'] = df['mutation'].map(codon_to_aa)
df['ref_aa'] = df['ref_codon'].map(codon_to_aa)
df = df.loc[df['Fragment']!='F1']

# Group by mutations
dgroup = dM.groupby(['Fragment','pcr','ref_codon','mut_codon','mutation','mut_pos','barcode']).agg(n_reads = ('seq','count')).reset_index()
dgroup['mut_code'] = dgroup.apply(lambda x: x['Fragment']+'-'+str(x['mut_codon'])+'-'+str(x['mutation']), axis=1)
dgroup['mut_codel'] = dgroup.apply(lambda x: x['Fragment']+'-'+str(x['mut_codon'])+'-'+str(x['mut_pos'])+'-'+str(x['mutation']), axis=1)
dgroup['mutation_aa'] = dgroup['mutation'].apply(lambda x: codon_to_aa[x])
dgroup['mut_code_aa'] = dgroup.apply(lambda x: x['Fragment']+'-'+str(x['mut_codon'])+'-'+str(x['mutation_aa']), axis=1)
dgroup['mut_codel_aa'] = dgroup.apply(lambda x: x['Fragment']+'-'+str(x['mut_codon'])+'-'+str(x['mut_pos'])+'-'+str(x['mutation_aa']), axis=1)

dout2 = dgroup.groupby(['Fragment','pcr','ref_codon','mut_codon','mutation','mut_code','mut_codel']).agg(
    n_barcodes_2 = ('barcode','count'),
    n_reads = ('n_reads','sum')).reset_index()

dout2.to_csv('02_combine_results_mutations_min' + str(barcode_cov) + '.tsv', sep='\t', header=True, index=False)

dc = dout2.copy()
dc['mutation_aa'] = dc['mutation'].map(codon_to_aa)
dc['ref_aa'] = dc['ref_codon'].map(codon_to_aa)
dc['pos'] = dc['mut_codel'].apply(lambda x: int(x.split('-')[2]))
dc['codon_pos'] = dc.apply(lambda x: abs(((x['mut_codon']*3) - 2) - x['pos']) if x['mut_codon']>0 else 0, axis=1)
dc['check'] = dc['mutation'].apply(lambda x: True if x[2] in ['G','T'] else False)
dc = dc[dc['check'] == True].reset_index(drop=True)
dc['nt'] = dc.apply(lambda x: x['mutation'][x['codon_pos']], axis=1)
dc = dc.loc[dc['Fragment']!='F1']
dc.head()

## Cumulative pcr analysis
# Generate 100 Random PCR subsets for resampling
pcrs = list(range(1,17))
i=10
S = []
for n in range(100):
    random.seed(i)
    sub = random.sample(pcrs,16)
    S.append(sub)
    i+=1

# Coverage (mutation in aa)
n=0
B = []

for pcr_set in S:
    print(n)

    R = []
    pcrs = []
    for npcr in range(len(pcr_set)):
        #print(pcr)
        pcr = pcr_set[npcr]
        pcrs.append(pcr)

        df_ = df[df['pcr'].isin(pcrs)]

        dgroup = df_.groupby(['Fragment','ref_aa','mut_codon','mutation_aa','barcode']).agg(n_reads = ('seq','count')).reset_index()
        dgroup['mut_code'] = dgroup.apply(lambda x: x['Fragment']+'-'+str(x['mut_codon'])+'-'+str(x['mutation_aa']), axis=1)

        dout = dgroup.groupby(['barcode','Fragment']).agg(n_muts = ('mut_code','nunique'),
                                                       n_ref = ('mut_code',calc_nunique_refs)).reset_index()
        dout['class'] = dout.apply(lambda x: "GOOD" if (x['n_muts']==1) & (x["n_ref"]==0) else "BAD", axis=1)

        df1 = dout.groupby(['Fragment'])['class'].value_counts().reset_index()
        df2 = dout.groupby(['Fragment']).agg(msum=('n_muts','count')).reset_index()
        df3 = pd.merge(df1, df2, on = ['Fragment'], how = 'left')
        df3['n_pcr'] = npcr
        df3['perc'] = df3.apply(lambda x: (x['count']/x['msum'])*100, axis=1)
        R.append(df3)
    dR = pd.concat(R, ignore_index=True)
    dR = dR[dR['class'] == 'GOOD'].reset_index(drop=True)
    dR.head()

    F = []
    for fragment in ['F1', 'F13', 'F43']:
        dc_frag = dc[dc['Fragment'] == fragment].reset_index(drop=True)
        C = []
        pcrs = []
        for npcr in range(len(pcr_set)):
            pcr = pcr_set[npcr]
            pcrs.append(pcr)
            
            dc_ = dc_frag[dc_frag['pcr'].isin(pcrs)].reset_index(drop=True)
            dcov = dc_.groupby(['Fragment','mut_codon']).agg(n_muts = ('mutation_aa','nunique')).reset_index()
            dcov['pcr'] = pcr
            dcov = dcov[dcov['mut_codon'] != 0].reset_index(drop=True)
            dres = dcov.groupby('Fragment').agg(mut_sum = ('n_muts','sum'), mut_avg = ('n_muts','mean')).reset_index()
            dres['n_pcr'] = npcr
            C.append(dres)
        dC = pd.concat(C, ignore_index=True)

        n_aa = dc_frag['mut_codon'].max()
        dC['tot_mut'] = n_aa * 21
        dC['perc_mutations'] = dC.apply(lambda x: (x['mut_sum']/x['tot_mut']) * 100, axis=1)
        F.append(dC)
    dF = pd.concat(F, ignore_index=True)

    dM = pd.merge(dR, dF, on = ['Fragment', 'n_pcr'], how = 'left')
    dM['pcr_set'] = n
    n+=1
    B.append(dM)
    
dB = pd.concat(B, ignore_index=True)
dB.to_csv("10_analyze_cumulative_barcodes_"+barcode_cov+"_aa.csv", sep=",", header=True, index=False)

# Barcode diversity (mutation in aa)
n=0
Z = []
for pcr_set in S:
    print(n)
    D = []
    pcrs = []
    for npcr in range(len(pcr_set)):
        #print(pcr)
        pcr = pcr_set[npcr]
        pcrs.append(pcr)
        df_ = df[df['pcr'].isin(pcrs)]

        dgroup = df_.groupby(['Fragment','ref_aa','mut_codon','mutation_aa','barcode']).agg(n_reads = ('seq','count')).reset_index()
        dgroup['mut_code'] = dgroup.apply(lambda x: x['Fragment']+'-'+str(x['mut_codon'])+'-'+str(x['mutation_aa']), axis=1)

        dout = dgroup.groupby(['Fragment','ref_aa','mut_codon','mutation_aa','mut_code']).agg(
            n_barcodes_2 = ('barcode','count'),
            reads = ('n_reads','sum')).reset_index()
        dout = dout[dout['mutation_aa'] != 'REF'].reset_index(drop=True)
        dout['barcode_4'] = dout['n_barcodes_2'].apply(lambda x: 1 if x>4 else 0)
        dout['barcode_9'] = dout['n_barcodes_2'].apply(lambda x: 1 if x>9 else 0)
        dg = dout.groupby(['Fragment']).agg(More_than_4 = ('barcode_4','sum'),
                                            More_than_9 = ('barcode_9','sum'),
                                            n = ('n_barcodes_2','count')).reset_index()
        dg['>4'] = dg.apply(lambda x: (x['More_than_4']/x['n'])*100, axis=1)
        dg['>9'] = dg.apply(lambda x: (x['More_than_9']/x['n'])*100, axis=1)

        dg['n_pcr'] = npcr
        D.append(dg)
    dD = pd.concat(D, ignore_index=True)
    dD['pcr_set'] = n
    n+=1
    Z.append(dD)
    
dZ = pd.concat(Z, ignore_index=True)
dZ.to_csv("10_analyze_cumulative_barcodes_diversity_"+barcode_cov+"_aa.csv", sep=",", header=True, index=False)   

# Combine results
dC = pd.merge(dB, dZ, on = ['Fragment','n_pcr','pcr_set'], how = 'left')
dC.loc[dC['Fragment']=='F1','nclones_pcr']=5000  
dC.loc[dC['Fragment']=='F13','nclones_pcr']=5000  
dC.loc[dC['Fragment']=='F43','nclones_pcr']=7500
dC['nclones']= dC['n_pcr']*dC['nclones_pcr']
dC.loc[dC['Fragment']=='F1','nclones_bp']= dC['nclones']/75
dC.loc[dC['Fragment']=='F13','nclones_bp']= dC['nclones']/75
dC.loc[dC['Fragment']=='F43','nclones_bp']= dC['nclones']/54
dC['nclones_bp'] = dC['nclones_bp'].round(0).astype('Int64')

# PLOT
sns.set_style("whitegrid")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (8, 5.33))

sns.lineplot(data = dC[dC['Fragment'] == 'F13'], x = 'nclones_bp', y = 'perc', ax=ax1, label="% informative barcodes", errorbar="ci")
sns.lineplot(data = dC[dC['Fragment'] == 'F13'], x = 'nclones_bp', y = 'perc_mutations', ax=ax1, label = "% mutation coverage", errorbar="ci")
sns.lineplot(data = dC[dC['Fragment'] == 'F13'], x = 'nclones_bp', y = '>4', ax=ax1, label = "% barcode diversity for mutations\nrepresented by >4 barcodes", errorbar="ci")
sns.lineplot(data = dC[dC['Fragment'] == 'F13'], x = 'nclones_bp', y = '>9', ax=ax1, label = "% barcode diversity for mutations\nrepresented by >9 barcodes", errorbar="ci")

sns.lineplot(data = dC[dC['Fragment'] == 'F43'], x = 'nclones_bp', y = 'perc', ax=ax2, label="% informative barcodes", errorbar="ci")
sns.lineplot(data = dC[dC['Fragment'] == 'F43'], x = 'nclones_bp', y = 'perc_mutations', ax=ax2, label = "% mutation coverage", errorbar="ci")
sns.lineplot(data = dC[dC['Fragment'] == 'F43'], x = 'nclones_bp', y = '>4', ax=ax2, label = "% barcode diversity for mutations\nrepresented by >4 barcodes", errorbar="ci")
sns.lineplot(data = dC[dC['Fragment'] == 'F43'], x = 'nclones_bp', y = '>9', ax=ax2, label = "% barcode diversity for mutations\nrepresented by >9 barcodes", errorbar="ci")

ax1.set_title("F13",fontsize=12)
ax2.set_title("F43",fontsize=12)

ax1.set_xlabel('Number of transformants per base pair (bp)', fontsize=12)
ax1.set_ylabel('Percentage (%)',fontsize=12)

ax2.set_xlabel('Number of transformants per base pair (bp)',fontsize=12)
ax2.set_ylabel('Percentage (%)',fontsize=12)

ax1.set_ylim(10,100)
ax2.set_ylim(10,100)
ax1.set_xlim(0,1000)
ax2.set_xlim(0,1000)

sns.move_legend(ax1, "upper left", bbox_to_anchor=(1, 1))
sns.move_legend(ax2, "upper left", bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig("10_analyze_cumulative_barcodes_vs_mutations_vs_diversity_100_random_nbclones_lineplot_1_aa.png", dpi=150)