#!/usr/bin/env python3
import csv
import sys
from Bio import SeqIO

# Parse command-line arguments
for i in range(len(sys.argv)):
    if sys.argv[i] == "--ref_seq": ref_seq = sys.argv[i+1]
    elif sys.argv[i] == "--consensus_seq": consensus_seq = sys.argv[i+1]
    elif sys.argv[i] == "--out_csv": out_csv = sys.argv[i+1]
    elif sys.argv[i] == "--sample_name": sample_name = sys.argv[i+1]
    elif sys.argv[i] == "--prot_name": prot_name = sys.argv[i+1]

aa_classes = {
    'A': 'Hydrophobic', 'V': 'Hydrophobic', 'L': 'Hydrophobic', 'M': 'Hydrophobic',
    'I': 'Hydrophobic', 'S': 'Polar non charged', 'T': 'Polar non charged',
    'N': 'Polar non charged', 'Q': 'Polar non charged', 'G': 'Special case',
    'C': 'Special case', 'P': 'Special case', 'U': 'Special case',
    'F': 'Aromatic Hydrophobic', 'Y': 'Aromatic Hydrophobic', 'W': 'Aromatic Hydrophobic',
    'K': 'Positively charged', 'R': 'Positively charged', 'H': 'Positively charged',
    'D': 'Positively charged', 'E': 'Positively charged', '*': 'Stop'
}

codon_code = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
}

def translate(seq):
    return ''.join(codon_code.get(seq[i:i+3], 'X') for i in range(0, len(seq)-2, 3))

def compare_sequences(ref_seq, cons_seq, sample_id, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(['SampleID', 'Protein', 'Mutation_type', 'Aa_change', 'Amino_Acid_Property_Change', 'Nt_mutation'])
        
        for codon_start in range(0, len(ref_seq)-2, 3):
            r_codon = ref_seq[codon_start:codon_start+3]
            c_codon = cons_seq[codon_start:codon_start+3]
            if r_codon != c_codon:
                mutations = [f"{r}{pos+codon_start+1}{c}" for pos, (r, c) in enumerate(zip(r_codon, c_codon)) if r != c]
                
                r_aa = codon_code.get(r_codon, '-')
                c_aa = codon_code.get(c_codon, '-')

                if '-' in r_codon:
                    if '-' in c_codon:
                        if r_codon.count('-') >= c_codon.count('-'):
                            # FRAMESHIFT
                            mutation_type = 'NON_SYNONYMOUS'
                            aa_desc = "Frameshift"
                            r_aa = "Framsehift"
                            c_aa = ""
                            
                    else:
                        # INSERTION
                        mutation_type = 'NON_SYNONYMOUS'
                        aa_desc = "Insertion"

                elif '-' in c_codon:
                    if c_codon.count('-') == 3:
                        # DELETION
                        mutation_type = 'NON_SYNONYMOUS'
                        aa_desc = "Deletion"
                    else:
                        # FRAMESHIFT
                        mutation_type = 'NON_SYNONYMOUS'
                        aa_desc = "Frameshift"
                        r_aa = "Framsehift"
                        c_aa = ""
                else:
                    if r_aa != c_aa:
                        mutation_type = 'NON_SYNONYMOUS'
                        aa_desc = f"Amino acid changed from {aa_classes.get(r_aa, 'unknown')} to {aa_classes.get(c_aa, 'unknown')}"
                    else:
                        mutation_type = 'SYNONYMOUS'
                        aa_desc = f"Amino acid did not change, it stayed {aa_classes.get(r_aa, 'unknown')}"

                writer.writerow([sample_id, prot_name, mutation_type, f"{r_aa}{(codon_start//3)+1}{c_aa}", aa_desc, ','.join(mutations)])

# Load sequences and run comparison
ref_record = list(SeqIO.parse(ref_seq, 'fasta'))[0].upper()
cons_record = list(SeqIO.parse(consensus_seq, 'fasta'))[0].upper()
compare_sequences(str(ref_record.seq), str(cons_record.seq), sample_name, out_csv)
