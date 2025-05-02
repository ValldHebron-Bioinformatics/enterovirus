import sys
import pandas as pd
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="Replace positions with coverage below 20 with 'N' in a FASTA sequence.")
parser.add_argument("--fasta", required=True, help="Path to the input FASTA file.")
parser.add_argument("--coverage", required=True, help="Path to the coverage CSV file.")
args = parser.parse_args()

coverage_file = pd.read_csv(args.coverage, sep=";")

fasta_records = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

# Process each sequence in the FASTA file
for ref, record in fasta_records.items():
    sequence = list(str(record.seq))
    ref_coverage = coverage_file[coverage_file["ref"] == ref]
    # Replace positions with coverage below 20 with "N"
    for _, row in ref_coverage.iterrows():
        if int(row["depth"]) < 20:
            sequence[int(row["pos"]) - 1] = "N"
    record.seq = Seq("".join(sequence))

# Write the modified sequences back to the FASTA file
with open(args.fasta, "w") as outfile:
    SeqIO.write(fasta_records.values(), outfile, "fasta")
