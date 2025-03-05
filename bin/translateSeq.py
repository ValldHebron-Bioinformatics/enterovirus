from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import sys

for i in range (len(sys.argv)):
    if sys.argv[i] == "--seq": seqFile = sys.argv[i+1]
    elif sys.argv[i] == "--fastaNucl": nuclFile = sys.argv[i+1]
    elif sys.argv[i] == "--fastaProt": protFile = sys.argv[i+1]

def translate_frames(nucleotide_seq):
    standard_table = CodonTable.unambiguous_dna_by_id[1]  # Standard genetic code
    stop_codons = standard_table.stop_codons
    start_codons = standard_table.start_codons
    
    longest_protein = ""
    longest_nucleotide_seq = ""
    
    for frame in range(3):
        seq = nucleotide_seq[frame:]
        protein = ""
        nucleotide_fragment = ""
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in stop_codons:
                if len(protein) > len(longest_protein):
                    longest_protein = protein
                    longest_nucleotide_seq = nucleotide_fragment
                protein = ""
                nucleotide_fragment = ""
            else:
                protein += Seq(codon).translate(table=standard_table)
                nucleotide_fragment += codon
        
        if len(protein) > len(longest_protein):
            longest_protein = protein
            longest_nucleotide_seq = nucleotide_fragment

    return longest_nucleotide_seq, longest_protein

def save_fasta(idSeq, sequence, filename):
    with open(filename, "a") as f:
        f.write(f">{idSeq}\n")
        f.write(sequence + "\n")

for seq_record in SeqIO.parse(seqFile, "fasta"):
    longest_nucleotide_seq, longest_protein = translate_frames(seq_record.seq)
    save_fasta(seq_record.id, str(longest_nucleotide_seq), nuclFile)
    save_fasta(seq_record.id, str(longest_protein), protFile)
