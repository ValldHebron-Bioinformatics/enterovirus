import csv
from Bio import SeqIO
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Calculate sequence length and percentage of Ns in a multi-FASTA file and append results to a CSV file.")
parser.add_argument("--fasta", required=True, help="Path to the input FASTA file.")
parser.add_argument("--csv", required=True, help="Path to the existing CSV file to append results.")
args = parser.parse_args()

# Open the CSV file in append mode
with open(args.csv, "a", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter=";")
    
    # Process each sequence in the FASTA file
    for record in SeqIO.parse(args.fasta, "fasta"):
        seq_length = len(record.seq)
        n_count = record.seq.upper().count("N")
        percentage_n = (n_count / seq_length) * 100
        # Write the results to the CSV file
        writer.writerow([f"{record.id}_length", seq_length])
        writer.writerow([f"{record.id}_percentageN", f"{percentage_n:.2f}%"])