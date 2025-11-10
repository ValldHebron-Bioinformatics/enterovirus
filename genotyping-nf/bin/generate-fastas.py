#!/usr/bin/env python3
import pandas as pd
import sys
from Bio import SeqIO
import re
import os
import csv

for i in range (len(sys.argv)):
    if sys.argv[i] == "--blast": blast_output = sys.argv[i+1]
    elif sys.argv[i] == "--scaffolds": scaffolds = sys.argv[i+1]
    elif sys.argv[i] == "--out-dir": out_dir = sys.argv[i+1]
    elif sys.argv[i] == "--protocol": protocol = sys.argv[i+1]
    elif sys.argv[i] == "--refs": refs = sys.argv[i+1]
    elif sys.argv[i] == "--input": inputval = sys.argv[i+1]


results = pd.DataFrame(columns=['match', 'seq','start','end','species','genotype'])
speciesTypedf = pd.read_csv(refs, sep=';')
speciesType = dict(zip(speciesTypedf['type'], speciesTypedf['species']))

# Import data
blastn = pd.read_csv(blast_output, sep = '\t', header=None) # Columns: idquery    idrefdb    score    evalue    x    x    x

# Create dictionary of scaffolds
sequences = {rec.id: rec.seq for rec in SeqIO.parse(scaffolds, "fasta")}

# Order rows by score (3), remove duplicates by column idquery (1): Get the best score for each NODE
blastn = blastn.sort_values([2], ascending=False)
# If duplicates in column idrefdb (2), means that there are several scaffolds for one organism (need to map to generate consensus)
blastn = blastn.drop_duplicates(subset=[0], keep='first')

if (protocol == "complete" or inputval == "fasta"):
    # Extract from original fasta the duplicates identified and map
    with open(out_dir+"/ev-match.fasta", 'w') as outfile:
        for _, row in blastn.iterrows():
            seq_id = str(row[0])  # Sequence ID
            segment = sequences[seq_id][row[4]-1:row[5]]
            if (row[4] <= row[5]) and (row[6] <= row[7]): # plus/plus
                outfile.write(f">{seq_id}\n{segment}\n")
                results.loc[len(results)] = [row[1], seq_id, row[6], row[7], speciesType[row[1].split('_')[1]], '']
            else: # plus/minus
                segmentrev=segment.replace("A","t").replace("C","g").replace("G","c").replace("T","a")
                segment=segmentrev[::-1].upper()
                outfile.write(f">{seq_id}\n{segment}\n")
                results.loc[len(results)] = [row[1], seq_id, row[7], row[6], speciesType[row[1].split('_')[1]], '']
elif protocol == "partial":
    # Get the one with higher cov
    blastn[['node','nodeval','len','lenval','cov','covval']] = blastn[0].str.split('_',expand=True)
    blastn.covval = pd.to_numeric(blastn.covval, errors='coerce')
    blastn = blastn.sort_values(['covval'], ascending=False)
    with open(out_dir+"/ev-match.fasta", 'w') as outfile:
        seq_id = str(blastn[0][0])
        segment = sequences[seq_id][blastn[4][0]-1:blastn[5][0]]
        if (blastn[4][0] <= blastn[5][0]) and (blastn[6][0] <= blastn[7][0]):
            outfile.write(f">{seq_id}\n{segment}\n")
            results.loc[len(results)] = [blastn[1][0].split('_')[1], str(blastn[0][0]), blastn[6][0], blastn[7][0], speciesType[blastn[1][0].split('_')[1]], ''] 
        # Bloque para cuando se usa como referencia de blastn la VP1. De momento usaremos todo EV, por lo que las posiciones no deben corregirse.
        #    if speciesType[blastn[1][0].split('_')[1]] == "Enterovirus alphacoxsackie": st = 2502
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus betacoxsackie": st = 2508
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus coxsackiepol": st = 2544
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus deconjuncti": st = 2418
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus alpharhino": st = 2441
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus betarhino": st = 2387
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus cerhino": st = 2469
        #    results.loc[len(results)] = [blastn[1][0].split('_')[1], str(blastn[0][0]), blastn[6][0]+st, blastn[7][0]+st, speciesType[blastn[1][0].split('_')[1]], ''] 
        else:
            segmentrev=segment.replace("A","t").replace("C","g").replace("G","c").replace("T","a")
            segment=segmentrev[::-1].upper()
            outfile.write(f">{seq_id}\n{segment}\n")
            results.loc[len(results)] = [blastn[1][0].split('_')[1], str(blastn[0][0]), blastn[7][0], blastn[6][0], speciesType[blastn[1][0].split('_')[1]], '']
        # Bloque para cuando se usa como referencia de blastn la VP1. De momento usaremos todo EV, por lo que las posiciones no deben corregirse.
        #    if speciesType[blastn[1][0].split('_')[1]] == "Enterovirus alphacoxsackie": st = 2502
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus betacoxsackie": st = 2508
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus coxsackiepol": st = 2544
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus deconjuncti": st = 2418
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus alpharhino": st = 2441
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus betarhino": st = 2387
        #    elif speciesType[blastn[1][0].split('_')[1]] == "Enterovirus cerhino": st = 2469
        #    results.loc[len(results)] = [blastn[1][0].split('_')[1], str(blastn[0][0]), blastn[7][0]+st, blastn[6][0]+st, speciesType[blastn[1][0].split('_')[1]], '']outfile.close()

results.to_csv(out_dir+"/species-assignment.csv", index=False)
