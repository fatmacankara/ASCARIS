from add_sequence import *
import pandas as pd
import numpy as np

def uniprotSequenceMatch(data):
    print('Retrieving UniProt sequences...\n')

    canonical_fasta = pd.DataFrame(columns=['uniprotID', 'uniprotSequence'])
    up_list = list(set(data['uniprotID'].to_list()))
    for i in range(len(up_list)):
        canonical_fasta.at[i, 'uniprotSequence'] = get_uniprot_seq(up_list[i])
        canonical_fasta.at[i, 'uniprotID'] = up_list[i]

    canonical_fasta = canonical_fasta.drop_duplicates()
    isoform_fasta = pd.DataFrame(columns=['uniprotID', 'isoformSequence'])
    iso_dict = []
    for i in range(len(up_list)):
        iso_dict.append(get_isoforms(up_list[i]))

    index = 0
    for i in iso_dict:
        for key, val in i.items():
            isoform_fasta.at[index, 'uniprotID'] = key
            isoform_fasta.at[index, 'isoformSequence'] = val
            index += 1
    isoform_fasta = isoform_fasta.drop_duplicates()

    for i in isoform_fasta.index:
        isoform_fasta.at[i, 'whichIsoform'] = isoform_fasta.at[i, 'uniprotID'][7:10].strip()
        isoform_fasta.at[i, 'uniprotID'] = isoform_fasta.at[i, 'uniprotID'][0:6]
    print('Sequence files created...\n')

    data = data.merge(canonical_fasta, on='uniprotID', how='left')
    data = data.replace({'': np.NaN, 'nan': np.NaN})
    data['whichIsoform'] = np.NaN
    data['wt_sequence_match'] = np.NaN
    not_match_in_uniprot = data[data.uniprotSequence.isna()]
    uniprot_matched = data[~data.uniprotSequence.isna()]

    return not_match_in_uniprot, uniprot_matched, canonical_fasta, isoform_fasta
