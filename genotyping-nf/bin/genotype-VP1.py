#!/usr/bin/env python3
import pandas as pd
import sys

for i in range (len(sys.argv)):
    if sys.argv[i] == "--sort-result": sortresult = sys.argv[i+1]
    elif sys.argv[i] == "--dir": dirpath= sys.argv[i+1]

assignmentdf = pd.read_csv(dirpath+'/species-assignment.csv', sep=',')
sortdf = pd.read_csv(sortresult, sep='\t')

sortdf["genotype_from_dataset"] = sortdf["dataset"].str.split("_").str[-1]
df_merged = assignmentdf.merge(sortdf[["seqName", "genotype_from_dataset"]], left_on="seq", right_on="seqName", how="left")
df_merged["genotype"] = df_merged["genotype_from_dataset"]

df_final = df_merged[assignmentdf.columns]

df_final.to_csv(dirpath+'/results/species-assignment.csv', index=False)