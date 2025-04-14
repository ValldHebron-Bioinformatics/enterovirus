#!/usr/bin/env python3
import os
import pandas as pd
import sys

for i in range(len(sys.argv)):
    if sys.argv[i] == '--muts_file':
        muts_file = sys.argv[i + 1]
    if sys.argv[i] == '--out_dir':
        out_dir = sys.argv[i + 1]
    elif sys.argv[i] == '--annotate':
        a_file = sys.argv[i + 1]
    elif sys.argv[i] == "--sample_name": sample_name = sys.argv[i+1]

df_muts_ALL = pd.read_csv(muts_file, sep=';')

annotate = pd.read_csv(a_file, sep=';')
annotated_df = pd.DataFrame()
annotated_df['Annotated'] = ''
annotated_df['Annotated_mutation'] = ''
more_than_one = []
one = []
for index in annotate.index:
    row = annotate.iloc[index].str.split('\t')
    dfr = pd.DataFrame(row)
    dfr = dfr.T
    dfr.columns = annotate.columns
    muts = dfr['mutation']
    if len(str(muts).split('+')) > 1:
        more_than_one.append(list(annotate.iloc[index].T))
    else:
        one.append(list(annotate.iloc[index].T))

if (not one) and (not more_than_one):
    exit()

if (one):
    one = pd.DataFrame(one)
    one.columns = annotate.columns
    for index in one.index:
        if (df_muts_ALL.Protein[df_muts_ALL.Protein == str(one.name[index])].shape[0] > 0): 
            if (df_muts_ALL[df_muts_ALL.Aa_change == str(one.mutation[index])].shape[0] > 0):
                row = df_muts_ALL[df_muts_ALL.Aa_change == str(one.mutation[index])].reset_index(drop=True)
                row.columns = df_muts_ALL.columns
                row['Annotated'] = ['yes']*row.shape[0]
                row['Annotated_mutation'] = [str(one.mutation[index])]*row.shape[0]
                if annotated_df.size > 0:
                    annotated_df = pd.concat([annotated_df, row])
                else:
                    annotated_df = row
            elif ('POI' in str(one.mutation[index])):
                aux_df = df_muts_ALL
                aux_df['Order'] = aux_df.Aa_change.str.extract(r'(\d+)', expand=False)
                aux_df['Order'] = pd.to_numeric(aux_df['Order'])

                aux_df2 = one
                aux_df2['Order'] = aux_df2.mutation.str.extract(r'(\d+)', expand=False)
                aux_df2['Order'] = pd.to_numeric(aux_df2['Order'])
                # print(aux_df[aux_df.Order == str(aux_df2.Order[index])])
                if (aux_df[aux_df.Order == aux_df2.Order[index]].shape[0] > 0):
                    row = aux_df[aux_df.Order == aux_df2.Order[index]].reset_index(drop=True)
                    row.columns = [
                        "SampleID", "Protein", "Mutation_type", "Aa_change", "Amino_Acid_Property_Change",
                        "Nt_mutation", "Order"
                    ]
                    row = row.drop(['Order'], axis=1)
                    row['Annotated'] = ['yes']*row.shape[0]
                    aux_row = aux_df[aux_df.Order == aux_df2.Order[index]]
                    value = aux_row.Aa_change.values[0]
                    row['Annotated_mutation'] = [value+" is in a POI"]*row.shape[0]
                    if annotated_df.size > 0:
                        annotated_df = pd.concat([annotated_df, row])
                    else:
                        annotated_df = row
                df_muts_ALL = df_muts_ALL.drop(['Order'], axis = 1)
                one = one.drop(['Order'], axis = 1)
            else: continue

if (more_than_one):
    more_than_one = pd.DataFrame(more_than_one)
    more_than_one.columns = annotate.columns
    for index in more_than_one.index:
        continue_next = False
        BY_GENE_MUTS=df_muts_ALL[df_muts_ALL.Protein == more_than_one.name[index]]
        for muts in more_than_one.mutation[index].split('+'):
            if (BY_GENE_MUTS[BY_GENE_MUTS.Aa_change == muts].shape[0] == 0):
                continue_next=True

        if (continue_next): continue

        for muts in more_than_one.mutation[index].split('+'):
            if (BY_GENE_MUTS[BY_GENE_MUTS.Aa_change == muts].shape[0] > 0):
                row = BY_GENE_MUTS[BY_GENE_MUTS.Aa_change == muts].reset_index(drop=True)
                row.columns = BY_GENE_MUTS.columns
                row['Annotated'] = ['yes']*row.shape[0]
                row['Annotated_mutation'] = more_than_one.mutation[index]*row.shape[0]
                if annotated_df.size > 0:
                    annotated_df = pd.concat([annotated_df, row])
                else:
                    annotated_df = row

annotated_df = annotated_df.reset_index(drop=True)
print(annotated_df)
# print('Annotated mutations found:')
# print(annotated_df)

df_muts_ALL.columns = [
    "SampleID", "Protein", "Mutation_type", "Aa_change",
    "Amino_Acid_Property_Change", "Nt_mutation"
]
df_muts_ALL['ColorCode'] = df_muts_ALL['Mutation_type']
df_muts_ALL.loc[df_muts_ALL['ColorCode'] == 'NON_SYNONYMOUS', 'ColorCode'] =  '#2243f5'
df_muts_ALL.loc[df_muts_ALL['ColorCode'] == 'SYNONYMOUS', 'ColorCode'] = '#b9c3fa'
df_muts_ALL['Type'] = df_muts_ALL['Mutation_type']
df_muts_ALL.loc[df_muts_ALL['Type'] == 'NON_SYNONYMOUS', 'Type'] =  'Non Synonymous Mutation'
df_muts_ALL.loc[df_muts_ALL['Type'] == 'SYNONYMOUS', 'Type'] = 'Synonymous Mutation'
df_muts_ALL['Order'] = df_muts_ALL.Aa_change.str.extract(r'(\d+)', expand=False)
df_muts_ALL['Order'] = pd.to_numeric(df_muts_ALL['Order'])
df_muts_ALL = df_muts_ALL.sort_values('Order')


df_muts_ALL['Annotated_mutation'] = df_muts_ALL.Aa_change
df_muts_ALL['Annotated'] = 'no'

annot_m = annotated_df
annot_m.columns = [
    "SampleID", "Protein", "Mutation_type", "Aa_change",
    "Amino_Acid_Property_Change", "Nt_mutation", "Annotated", "Annotated_mutation"
]
annot_m['Order'] = annot_m.Aa_change.str.extract(r'(\d+)', expand=False)
annot_m['Order'] = pd.to_numeric(annot_m['Order'])
annot_m = annot_m.sort_values('Order')

for i in annot_m.Aa_change:
    df_muts_ALL.loc[df_muts_ALL.Aa_change == i, 'Annotated'] = 'yes'
    df_muts_ALL.loc[df_muts_ALL.Aa_change == i, 'Annotated_mutation'] = annot_m.loc[annot_m.Aa_change == i, 'Annotated_mutation'].values[0]


# print(df_muts_ALL)
df_muts_ALL.to_csv(out_dir + '/' + f'Annotated_mutations_{sample_name}.csv', sep=';', index=False)