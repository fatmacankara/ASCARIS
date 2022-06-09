import ssl
import requests as r
from decimal import *
import numpy as np
def add_annotations(dataframe):
    print('Downloading UniProt sequence annotations...\n')
    ssl._create_default_https_context = ssl._create_unverified_context

    original_annot_name = ['DISULFID', 'INIT_MET', 'INTRAMEM', 'VARIANT', 'DNA_BIND', 'ACT_SITE', 'NP_BIND', 'LIPID',
                           'SITE',
                           'TRANSMEM', 'CROSSLNK', 'MUTAGEN', 'STRAND', 'HELIX', 'TURN', 'METAL', 'REPEAT', 'TOPO_DOM',
                           'CA_BIND', 'BINDING', 'REGION', 'SIGNAL', 'MOD_RES', 'ZN_FING', 'MOTIF', 'COILED', 'PEPTIDE',
                           'TRANSIT', 'CARBOHYD', 'PROPEP']
    annotation_list = ['disulfide', 'intMet', 'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                       'nucleotideBinding', 'lipidation', 'site', 'transmembrane', 'crosslink', 'mutagenesis', 'strand',
                       'helix', 'turn', 'metalBinding', 'repeat', 'topologicalDomain', 'caBinding', 'bindingSite',
                       'region',
                       'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil', 'peptide',
                       'transitPeptide', 'glycosylation', 'propeptide']

    dataframe = dataframe.reset_index().drop(['index'], axis=1)

    for annot in original_annot_name:
        dataframe[annot] = ''

    for protein in list(set(dataframe.uniprotID.to_list())):
        print('Downloading annotations for ' + protein)
        #uniprot_entry = urllib.request.urlopen("http://www.uniprot.org/uniprot/" + protein + ".txt").read()
        #uniprot_entry = uniprot_entry.decode('utf-8').split('\n')
        uniprot_entry = r.get("http://www.uniprot.org/uniprot/" + protein + ".txt")
        uniprot_entry = uniprot_entry.text.split('\n')

        annot_for_protein = []
        for annotation in original_annot_name:
            for line in uniprot_entry:
                if annotation.strip() in line and line.startswith(
                        'FT') and 'evidence' not in line and 'ECO' not in line and 'note' not in line:
                    annot_for_protein.append(list(filter(None, line.split(' ')))[1:])
        for select in annot_for_protein:
            if select[0] not in dataframe.columns:
                dataframe.loc[dataframe.uniprotID == protein, select[0]] = str((select[1] + '; '))
            else:
                dataframe.loc[dataframe.uniprotID == protein, select[0]] += str((select[1] + '; '))

    for i in range(len(original_annot_name)):
        dataframe = dataframe.rename(columns={original_annot_name[i]: annotation_list[i]})

    # Fix annotation positions

    print()
    print('Processing positions...\n')
    for i in dataframe.index:
        for annot in dataframe.columns[-30:]:
            if annot != 'disulfide':
                if dataframe.at[i, annot] != 'nan':
                    dataframe.at[i, annot] = ([x for x in [k.strip() for k in dataframe.at[i, annot].split(';')] if x])
                    if '..' not in str(dataframe.at[i, annot]):
                        pass
                    elif '..' in str(dataframe.at[i, annot]):
                        dataframe.at[i, annot] = str(dataframe.at[i, annot]).replace('..', '-')
            else:
                disulfide_annot = []
                if dataframe.at[i, annot] != 'nan':
                    dataframe.at[i, annot]=  dataframe.at[i, annot].split(';')
                    dataframe.at[i, annot] = [i.split('..') for i in dataframe.at[i, annot]]
                    dataframe.at[i, annot] =[e for v in  dataframe.at[i, annot] for e in v]
                    dataframe.at[i, annot] = [i for i in dataframe.at[i, annot] if i != ' ']

    # Add binary annotations
    print('Adding binary annotations...\n')
    dataframe = dataframe.astype('str')
    for i in dataframe.index:
        for k in annotation_list:  # get the positions of each attribute as a list
            txt = k + 'Binary'
            dataframe.at[i, txt] = Decimal('nan')
            try:
                for positions in dataframe.at[i, k].split(','):
                    position = positions.strip('[').strip(']').replace("'", "")
                    if position != 'nan' and position != '' and '-' not in position and int(
                            dataframe.at[i, 'pos']) == int(position):
                        dataframe.at[i, txt] = '1'
                        break
                    elif position != 'nan' and position != '' and '-' not in position and int(
                            dataframe.at[i, 'pos']) != int(position):
                        dataframe.at[i, txt] = '0'
                    elif position != 'nan' and position != '' and '-' in position:
                        if int(position.split('-')[0]) < int(dataframe.at[i, 'pos']) < int(position.split('-')[1]):
                            # print(position.split('-')[0], position.split('-')[1])
                            dataframe.at[i, txt] = '1'
                            break
                        else:
                            dataframe.at[i, txt] = '0'
            except:
                ValueError

    # Final corrections

    dataframe = dataframe.replace({'[\'?\']': 'nan'})
    dataframe = dataframe.replace({'[]': 'nan'})
    return dataframe

