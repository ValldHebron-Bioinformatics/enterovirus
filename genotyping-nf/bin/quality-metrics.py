import csv
import pandas as pd
from Bio import SeqIO
import argparse

# Argument parser
parser = argparse.ArgumentParser(description = "Calculate sequence length and percentage of Ns in a multi-FASTA file and append results to a CSV file.")
parser.add_argument("--fasta", required = True, help = "Path to the input FASTA file.")
parser.add_argument("--csv", required = True, help = "Path to the existing CSV file to append results.")
parser.add_argument("--coverage", required = True, help = "Path to the coverage CSV file.")
args = parser.parse_args()

# Open the CSV file in append mode
with open(args.csv, "a", newline = "") as csvfile:
    writer = csv.writer(csvfile, delimiter = ";")
    cov_file = pd.read_csv(args.coverage, sep = ";")
    # Process each sequence in the FASTA file
    for record in SeqIO.parse(args.fasta, "fasta"):
        cov_record = cov_file[cov_file.ref == record.id]
        depth_coverage = cov_record.depth.quantile([0.5]).iloc[0]
        seq_length = len(record.seq)
        n_count = record.seq.upper().count("N")
        percentage_n = (n_count / seq_length) * 100
        # Write the results to the CSV file
        writer.writerow([f"{record.id}_length", seq_length])
        writer.writerow([f"{record.id}_percentageN", f"{percentage_n:.2f}%"])
        writer.writerow([f"{record.id}_median_depth_coverage", f"{depth_coverage:.2f}x"])
