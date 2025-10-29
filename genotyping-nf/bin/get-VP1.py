#!/usr/bin/env python3
import pandas as pd
import sys
from Bio import SeqIO
import csv

for i in range (len(sys.argv)):
    if sys.argv[i] == "--diamond": diamond = sys.argv[i+1]
    elif sys.argv[i] == "--dir": dirpath= sys.argv[i+1]
    elif sys.argv[i] == "--refs": refs = sys.argv[i+1]
    elif sys.argv[i] == "--pwd": pwd = sys.argv[i+1]


results = pd.DataFrame(columns=['SAMPLE_DIR', 'prot','VP1consensus','EVreference', "genotype"])
speciesTypedf = pd.read_csv(refs, sep=';')
speciesType = dict(zip(speciesTypedf['type'], speciesTypedf['species']))
blastndf = pd.read_csv(dirpath+'/species-assignment.csv', sep=',')

# Import data
diamonddf = pd.read_csv(diamond, sep = '\t', header=None)
diamonddf = diamonddf.sort_values([2], ascending=False)
diamonddf = diamonddf.drop_duplicates(subset=[0], keep='first')

# Load protein sequences
protein_sequences = {rec.id: rec.seq for rec in SeqIO.parse(dirpath+"/ev-match_cds-aa.fasta", "fasta")}
# Load nucleotide sequences
nucleotide_sequences = {rec.id: rec.seq for rec in SeqIO.parse(dirpath+"/ev-match_cds-nucl.fasta", "fasta")}

prot_output = open(dirpath+"/VP1_prot.fasta", "w")
nucl_output = open(dirpath+"/VP1_nucl.fasta", "w")

for _, row in diamonddf.iterrows():
    seq_id = row[0]  # Sequence ID
    start = row[4] # Adjusted start position
    end = row[5]  # Adjusted end position

    genotype = row[1].split('_')[1]
    blastndf.loc[blastndf['seq'] == seq_id, 'genotype'] = genotype

    # Extract protein sequence
    if seq_id in protein_sequences:
        prot_segment = protein_sequences[seq_id][start-1:end]
        prot_output.write(f">{seq_id}\n{prot_segment}\n")
    
    # Extract nucleotide sequence
    if seq_id in nucleotide_sequences:
        nucl_segment = nucleotide_sequences[seq_id][start*3-3:end*3]  # Convert amino acid positions to nucleotide
        nucl_output.write(f">{seq_id}\n{nucl_segment}\n")
    
    if speciesType[row[1].split('_')[1]] == "Enterovirus alphacoxsackie": sp = "EV-A_reference-VP1_nucleotide.fasta"
    elif speciesType[row[1].split('_')[1]] == "Enterovirus betacoxsackie": sp = "EV-B_reference-VP1_nucleotide.fasta"
    elif speciesType[row[1].split('_')[1]] == "Enterovirus coxsackiepol": sp = "EV-C_reference-VP1_nucleotide.fasta"
    elif speciesType[row[1].split('_')[1]] == "Enterovirus deconjuncti": sp = "EV-D_reference-VP1_nucleotide.fasta"
    results.loc[len(results)] = [dirpath, 'VP1', dirpath+"/"+seq_id+".fasta", pwd+"/"+sp, genotype]

# Close output files
prot_output.close()
nucl_output.close()

blastndf.to_csv(dirpath+"/species-assignment.csv", index=False)
results.to_csv(dirpath+"/sample-mutations.csv", index=False)
