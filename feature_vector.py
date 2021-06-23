# path_to_output_files  (output)
# path_to_domains (input)  bunlar verildi, ama input al ve koddan sil.
# swiss_model_path (input)
# path_to_interfaces
# path to fisher


# STEP 0 : IMPORT NECESSARY MODULES AND LIBRARIES
from timeit import default_timer as timer
import xml.etree.ElementTree as ET
from collections import Counter
from bs4 import BeautifulSoup
from io import StringIO
from decimal import *
import pandas as pd
import requests as r
import os.path as op
import subprocess
import ssbio.utils
import warnings
import sys
import pathlib
import os, glob
import math
import ssbio
import ssl
from Bio.Align import substitution_matrices
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBList
from Bio import Align
from Bio import SeqIO
from Bio.PDB import *

warnings.filterwarnings("ignore")
start = timer()

# FUNCTIONS

def compositionValues(aa1, aa2):
    compositionValues = {'S': 1.42, 'R': 0.65, 'L': 0, 'P': 0.39, 'T': 0.71, 'A': 0, 'V': 0, 'G': 0.74,
                         'I': 0, 'F': 0, 'Y': 0.20, 'C': 2.75, 'H': 0.58, 'Q': 0.89, 'N': 1.33, 'K': 0.33,
                         'D': 1.38, 'E': 0.92, 'M': 0, 'W': 0.13}
    dif = round((compositionValues[aa1] - compositionValues[aa2]), 2)
    return (dif)


def polarityValues(aa1, aa2):
    polarityValues = {'S': 9.2, 'R': 10.5, 'L': 4.9, 'P': 8.0, 'T': 8.6, 'A': 8.1, 'V': 5.9, 'G': 9.0,
                      'I': 5.2, 'F': 5.2, 'Y': 6.2, 'C': 5.5, 'H': 10.4, 'Q': 10.5, 'N': 11.6, 'K': 11.3,
                      'D': 13.0, 'E': 12.3, 'M': 5.7, 'W': 5.4}
    dif = round((polarityValues[aa1] - polarityValues[aa2]), 2)
    return (dif)


def volumeValues(aa1, aa2):
    volumeValues = {'S': 32, 'R': 124, 'L': 111, 'P': 32.5, 'T': 61, 'A': 31, 'V': 84, 'G': 3,
                    'I': 111, 'F': 132, 'Y': 136, 'C': 55, 'H': 96, 'Q': 85, 'N': 56, 'K': 119,
                    'D': 54, 'E': 83, 'M': 105, 'W': 170}
    dif = round((volumeValues[aa1] - volumeValues[aa2]), 2)
    return (dif)


def add_physicochemical(df):
    grantham_dict = {
        ('A', 'A'): '0',
        ('A', 'C'): '195',
        ('A', 'D'): '126',
        ('A', 'E'): '107',
        ('A', 'F'): '113',
        ('A', 'G'): '60',
        ('A', 'H'): '86',
        ('A', 'I'): '94',
        ('A', 'K'): '106',
        ('A', 'L'): '96',
        ('A', 'M'): '84',
        ('A', 'N'): '111',
        ('A', 'P'): '27',
        ('A', 'Q'): '91',
        ('A', 'R'): '112',
        ('A', 'S'): '99',
        ('A', 'T'): '58',
        ('A', 'V'): '64',
        ('A', 'W'): '148',
        ('A', 'Y'): '112',
        ('C', 'A'): '195',
        ('C', 'C'): '0',
        ('C', 'D'): '154',
        ('C', 'E'): '170',
        ('C', 'F'): '205',
        ('C', 'G'): '159',
        ('C', 'H'): '174',
        ('C', 'I'): '198',
        ('C', 'K'): '202',
        ('C', 'L'): '198',
        ('C', 'M'): '196',
        ('C', 'N'): '139',
        ('C', 'P'): '169',
        ('C', 'Q'): '154',
        ('C', 'R'): '180',
        ('C', 'S'): '112',
        ('C', 'T'): '149',
        ('C', 'V'): '192',
        ('C', 'W'): '215',
        ('C', 'Y'): '194',
        ('D', 'A'): '126',
        ('D', 'C'): '154',
        ('D', 'D'): '0',
        ('D', 'E'): '45',
        ('D', 'F'): '177',
        ('D', 'G'): '94',
        ('D', 'H'): '81',
        ('D', 'I'): '168',
        ('D', 'K'): '101',
        ('D', 'L'): '172',
        ('D', 'M'): '160',
        ('D', 'N'): '23',
        ('D', 'P'): '108',
        ('D', 'Q'): '61',
        ('D', 'R'): '96',
        ('D', 'S'): '65',
        ('D', 'T'): '85',
        ('D', 'V'): '152',
        ('D', 'W'): '181',
        ('D', 'Y'): '160',
        ('E', 'A'): '107',
        ('E', 'C'): '170',
        ('E', 'D'): '45',
        ('E', 'E'): '0',
        ('E', 'F'): '140',
        ('E', 'G'): '98',
        ('E', 'H'): '40',
        ('E', 'I'): '134',
        ('E', 'K'): '56',
        ('E', 'L'): '138',
        ('E', 'M'): '126',
        ('E', 'N'): '42',
        ('E', 'P'): '93',
        ('E', 'Q'): '29',
        ('E', 'R'): '54',
        ('E', 'S'): '80',
        ('E', 'T'): '65',
        ('E', 'V'): '121',
        ('E', 'W'): '152',
        ('E', 'Y'): '122',
        ('F', 'A'): '113',
        ('F', 'C'): '205',
        ('F', 'D'): '177',
        ('F', 'E'): '140',
        ('F', 'F'): '0',
        ('F', 'G'): '153',
        ('F', 'H'): '100',
        ('F', 'I'): '21',
        ('F', 'K'): '102',
        ('F', 'L'): '22',
        ('F', 'M'): '28',
        ('F', 'N'): '158',
        ('F', 'P'): '114',
        ('F', 'Q'): '116',
        ('F', 'R'): '97',
        ('F', 'S'): '155',
        ('F', 'T'): '103',
        ('F', 'V'): '50',
        ('F', 'W'): '40',
        ('F', 'Y'): '22',
        ('G', 'A'): '60',
        ('G', 'C'): '159',
        ('G', 'D'): '94',
        ('G', 'E'): '98',
        ('G', 'F'): '153',
        ('G', 'G'): '0',
        ('G', 'H'): '98',
        ('G', 'I'): '135',
        ('G', 'K'): '127',
        ('G', 'L'): '138',
        ('G', 'M'): '127',
        ('G', 'N'): '80',
        ('G', 'P'): '42',
        ('G', 'Q'): '87',
        ('G', 'R'): '125',
        ('G', 'S'): '56',
        ('G', 'T'): '59',
        ('G', 'V'): '109',
        ('G', 'W'): '184',
        ('G', 'Y'): '147',
        ('H', 'A'): '86',
        ('H', 'C'): '174',
        ('H', 'D'): '81',
        ('H', 'E'): '40',
        ('H', 'F'): '100',
        ('H', 'G'): '98',
        ('H', 'H'): '0',
        ('H', 'I'): '94',
        ('H', 'K'): '32',
        ('H', 'L'): '99',
        ('H', 'M'): '87',
        ('H', 'N'): '68',
        ('H', 'P'): '77',
        ('H', 'Q'): '24',
        ('H', 'R'): '29',
        ('H', 'S'): '89',
        ('H', 'T'): '47',
        ('H', 'V'): '84',
        ('H', 'W'): '115',
        ('H', 'Y'): '83',
        ('I', 'A'): '94',
        ('I', 'C'): '198',
        ('I', 'D'): '168',
        ('I', 'E'): '134',
        ('I', 'F'): '21',
        ('I', 'G'): '135',
        ('I', 'H'): '94',
        ('I', 'I'): '0',
        ('I', 'K'): '102',
        ('I', 'L'): '5',
        ('I', 'M'): '10',
        ('I', 'N'): '149',
        ('I', 'P'): '95',
        ('I', 'Q'): '109',
        ('I', 'R'): '97',
        ('I', 'S'): '142',
        ('I', 'T'): '89',
        ('I', 'V'): '29',
        ('I', 'W'): '61',
        ('I', 'Y'): '33',
        ('K', 'A'): '106',
        ('K', 'C'): '202',
        ('K', 'D'): '101',
        ('K', 'E'): '56',
        ('K', 'F'): '102',
        ('K', 'G'): '127',
        ('K', 'H'): '32',
        ('K', 'I'): '102',
        ('K', 'K'): '0',
        ('K', 'L'): '107',
        ('K', 'M'): '95',
        ('K', 'N'): '94',
        ('K', 'P'): '103',
        ('K', 'Q'): '53',
        ('K', 'R'): '26',
        ('K', 'S'): '121',
        ('K', 'T'): '78',
        ('K', 'V'): '97',
        ('K', 'W'): '110',
        ('K', 'Y'): '85',
        ('L', 'A'): '96',
        ('L', 'C'): '198',
        ('L', 'D'): '172',
        ('L', 'E'): '138',
        ('L', 'F'): '22',
        ('L', 'G'): '138',
        ('L', 'H'): '99',
        ('L', 'I'): '5',
        ('L', 'K'): '107',
        ('L', 'L'): '0',
        ('L', 'M'): '15',
        ('L', 'N'): '153',
        ('L', 'P'): '98',
        ('L', 'Q'): '113',
        ('L', 'R'): '102',
        ('L', 'S'): '145',
        ('L', 'T'): '92',
        ('L', 'V'): '32',
        ('L', 'W'): '61',
        ('L', 'Y'): '36',
        ('M', 'A'): '84',
        ('M', 'C'): '196',
        ('M', 'D'): '160',
        ('M', 'E'): '126',
        ('M', 'F'): '28',
        ('M', 'G'): '127',
        ('M', 'H'): '87',
        ('M', 'I'): '10',
        ('M', 'K'): '95',
        ('M', 'L'): '15',
        ('M', 'M'): '0',
        ('M', 'N'): '142',
        ('M', 'P'): '87',
        ('M', 'Q'): '101',
        ('M', 'R'): '91',
        ('M', 'S'): '135',
        ('M', 'T'): '81',
        ('M', 'V'): '21',
        ('M', 'W'): '67',
        ('M', 'Y'): '36',
        ('N', 'A'): '111',
        ('N', 'C'): '139',
        ('N', 'D'): '23',
        ('N', 'E'): '42',
        ('N', 'F'): '158',
        ('N', 'G'): '80',
        ('N', 'H'): '68',
        ('N', 'I'): '149',
        ('N', 'K'): '94',
        ('N', 'L'): '153',
        ('N', 'M'): '142',
        ('N', 'N'): '0',
        ('N', 'P'): '91',
        ('N', 'Q'): '46',
        ('N', 'R'): '86',
        ('N', 'S'): '46',
        ('N', 'T'): '65',
        ('N', 'V'): '133',
        ('N', 'W'): '174',
        ('N', 'Y'): '143',
        ('P', 'A'): '27',
        ('P', 'C'): '169',
        ('P', 'D'): '108',
        ('P', 'E'): '93',
        ('P', 'F'): '114',
        ('P', 'G'): '42',
        ('P', 'H'): '77',
        ('P', 'I'): '95',
        ('P', 'K'): '103',
        ('P', 'L'): '98',
        ('P', 'M'): '87',
        ('P', 'N'): '91',
        ('P', 'P'): '0',
        ('P', 'Q'): '76',
        ('P', 'R'): '103',
        ('P', 'S'): '74',
        ('P', 'T'): '38',
        ('P', 'V'): '68',
        ('P', 'W'): '147',
        ('P', 'Y'): '110',
        ('Q', 'A'): '91',
        ('Q', 'C'): '154',
        ('Q', 'D'): '61',
        ('Q', 'E'): '29',
        ('Q', 'F'): '116',
        ('Q', 'G'): '87',
        ('Q', 'H'): '24',
        ('Q', 'I'): '109',
        ('Q', 'K'): '53',
        ('Q', 'L'): '113',
        ('Q', 'M'): '101',
        ('Q', 'N'): '46',
        ('Q', 'P'): '76',
        ('Q', 'Q'): '0',
        ('Q', 'R'): '43',
        ('Q', 'S'): '68',
        ('Q', 'T'): '42',
        ('Q', 'V'): '96',
        ('Q', 'W'): '130',
        ('Q', 'Y'): '99',
        ('R', 'A'): '112',
        ('R', 'C'): '180',
        ('R', 'D'): '96',
        ('R', 'E'): '54',
        ('R', 'F'): '97',
        ('R', 'G'): '125',
        ('R', 'H'): '29',
        ('R', 'I'): '97',
        ('R', 'K'): '26',
        ('R', 'L'): '102',
        ('R', 'M'): '91',
        ('R', 'N'): '86',
        ('R', 'P'): '103',
        ('R', 'Q'): '43',
        ('R', 'R'): '0',
        ('R', 'S'): '110',
        ('R', 'T'): '71',
        ('R', 'V'): '96',
        ('R', 'W'): '101',
        ('R', 'Y'): '77',
        ('S', 'A'): '99',
        ('S', 'C'): '112',
        ('S', 'D'): '65',
        ('S', 'E'): '80',
        ('S', 'F'): '155',
        ('S', 'G'): '56',
        ('S', 'H'): '89',
        ('S', 'I'): '142',
        ('S', 'K'): '121',
        ('S', 'L'): '145',
        ('S', 'M'): '135',
        ('S', 'N'): '46',
        ('S', 'P'): '74',
        ('S', 'Q'): '68',
        ('S', 'R'): '110',
        ('S', 'S'): '0',
        ('S', 'T'): '58',
        ('S', 'V'): '124',
        ('S', 'W'): '177',
        ('S', 'Y'): '144',
        ('T', 'A'): '58',
        ('T', 'C'): '149',
        ('T', 'D'): '85',
        ('T', 'E'): '65',
        ('T', 'F'): '103',
        ('T', 'G'): '59',
        ('T', 'H'): '47',
        ('T', 'I'): '89',
        ('T', 'K'): '78',
        ('T', 'L'): '92',
        ('T', 'M'): '81',
        ('T', 'N'): '65',
        ('T', 'P'): '38',
        ('T', 'Q'): '42',
        ('T', 'R'): '71',
        ('T', 'S'): '58',
        ('T', 'T'): '0',
        ('T', 'V'): '69',
        ('T', 'W'): '128',
        ('T', 'Y'): '92',
        ('V', 'A'): '64',
        ('V', 'C'): '192',
        ('V', 'D'): '152',
        ('V', 'E'): '121',
        ('V', 'F'): '50',
        ('V', 'G'): '109',
        ('V', 'H'): '84',
        ('V', 'I'): '29',
        ('V', 'K'): '97',
        ('V', 'L'): '32',
        ('V', 'M'): '21',
        ('V', 'N'): '133',
        ('V', 'P'): '68',
        ('V', 'Q'): '96',
        ('V', 'R'): '96',
        ('V', 'S'): '124',
        ('V', 'T'): '69',
        ('V', 'V'): '0',
        ('V', 'W'): '88',
        ('V', 'Y'): '55',
        ('W', 'A'): '148',
        ('W', 'C'): '215',
        ('W', 'D'): '181',
        ('W', 'E'): '152',
        ('W', 'F'): '40',
        ('W', 'G'): '184',
        ('W', 'H'): '115',
        ('W', 'I'): '61',
        ('W', 'K'): '110',
        ('W', 'L'): '61',
        ('W', 'M'): '67',
        ('W', 'N'): '174',
        ('W', 'P'): '147',
        ('W', 'Q'): '130',
        ('W', 'R'): '101',
        ('W', 'S'): '177',
        ('W', 'T'): '128',
        ('W', 'V'): '88',
        ('W', 'W'): '0',
        ('W', 'Y'): '37',
        ('Y', 'A'): '112',
        ('Y', 'C'): '194',
        ('Y', 'D'): '160',
        ('Y', 'E'): '122',
        ('Y', 'F'): '22',
        ('Y', 'G'): '147',
        ('Y', 'H'): '83',
        ('Y', 'I'): '33',
        ('Y', 'K'): '85',
        ('Y', 'L'): '36',
        ('Y', 'M'): '36',
        ('Y', 'N'): '143',
        ('Y', 'P'): '110',
        ('Y', 'Q'): '99',
        ('Y', 'R'): '77',
        ('Y', 'S'): '144',
        ('Y', 'T'): '92',
        ('Y', 'V'): '55',
        ('Y', 'W'): '37',
        ('Y', 'Y'): '0'
    }
    for i in df.index:
        try:
            df.at[i, 'composition'] = compositionValues(df.at[i, 'wt'], df.at[i, 'mut'])
            df.at[i, 'polarity'] = polarityValues(df.at[i, 'wt'], df.at[i, 'mut'])
            df.at[i, 'volume'] = volumeValues(df.at[i, 'wt'], df.at[i, 'mut'])
            df.at[i, 'granthamScore'] = grantham_dict[df.at[i, 'wt'], df.at[i, 'mut']]
        except:
            KeyError
            df.at[i, 'composition'] = 'nan'
            df.at[i, 'polarity'] = 'nan'
            df.at[i, 'volume'] = 'nan'
            df.at[i, 'granthamScore'] = 'nan'
    return df


def add_domains(data, path_to_domains):
    domains = pd.read_csv(path_to_domains, delimiter=' ')
    data = data.merge(domains, right_on='proteinID', left_on='uniprotID', how='left')
    data = data.drop(['proteinID'], axis=1)
    # Label each data point as range or notRange based on the relative distance of mutation and domain boundaries.
    data = data.astype('str')
    data.domStart = data.domStart.astype('float')
    data.domEnd = data.domEnd.astype('float')

    for i in data.index:
        if data.at[i, 'domain'] != 'nan':
            if int(data.at[i, 'domStart']) <= int(data.at[i, 'pos']) <= int(data.at[i, 'domEnd']):
                data.at[i, 'distance'] = 0
            else:
                distance = min(abs(int(data.at[i, 'domStart']) - int(data.at[i, 'pos'])),
                               abs(int(data.at[i, 'domEnd']) - int(data.at[i, 'pos'])))
                data.at[i, 'distance'] = int(distance)
        else:
            data.at[i, 'distance'] = 'nan'

    data = data.sort_values(by=['datapoint', 'distance']).reset_index(drop=True)  # Distances will be sorted.

    # Keep the one with the least distance. But we may have more than one range domains for a datapoint if distance = 0.
    # For this reason first we need to separate range ones so that when we take the first occurance to get the closest one
    # for non range ones, other distance=0 ones wont disappear.

    data_range = data[data.distance == 0]
    data_out_range = data[data.distance != 0]

    # For the range ones, find the most occurance

    dom = []
    for i in data_range.index:
        dom.append(data_range.at[i, 'domain'])

    domainCount = Counter(dom)  # Occurance of domains.

    # For out of range ones, take the closest distance.
    data_out_range = data_out_range.drop_duplicates(['datapoint'], keep='first')  # Already sorted above.
    domain_counts = pd.DataFrame(domainCount.items(), columns=['domain', 'count'])
    data_range_counts = data_range.merge(domain_counts, on='domain')
    data_range_counts = data_range_counts.sort_values(['datapoint', 'count'])
    data_range_counts = data_range_counts.drop_duplicates(['datapoint'], keep='last')  # Take with the higher count.
    data_range_counts = data_range_counts.drop(['count'], axis=1)

    # Merge them back together

    frames = [data_range_counts, data_out_range]
    data = pd.concat(frames, sort=False)  # Here when you concat two data frames, we might have range and not range with
    # min distance for the same data point. Delete the one coming from notRange one.
    data = data.sort_values(['datapoint', 'distance']).reset_index(drop=True)
    data = data.drop_duplicates(['datapoint'], keep='first')

    data = data[['uniprotID', 'wt', 'pos', 'mut', 'composition', 'polarity',
                 'volume', 'domain', 'domStart', 'domEnd', 'distance', 'datapoint']]
    data = data.astype(str)
    return data


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
            if dataframe.at[i, annot] != 'nan':
                dataframe.at[i, annot] = ([x for x in [k.strip() for k in dataframe.at[i, annot].split(';')] if x])
                if '..' not in str(dataframe.at[i, annot]):
                    pass
                elif '..' in str(dataframe.at[i, annot]):
                    dataframe.at[i, annot] = str(dataframe.at[i, annot]).replace('..', '-')

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


def get_uniprot_seq(protein_id):
    print('Fetching for ID: ', protein_id)
    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl + protein_id + ".fasta"
    response = r.post(currentUrl)
    cData = ''.join(response.text)
    Seq = StringIO(cData)
    pSeq = list(SeqIO.parse(Seq, 'fasta'))
    try:
        return str(pSeq[0].seq)
    except:
        IndexError
        return str('')


def get_isoforms(protein_id):
    print('Fetching for ID: ', protein_id)
    try:
        # a dictionary storing the sequence of your isoforms, key: accesion number, value: sequence
        isoforms = dict()
        # make a call to EBI API
        req = r.get('https://www.ebi.ac.uk/proteins/api/proteins/{}/isoforms'.format(protein_id))
        # parse the returned XML
        uniprot = ET.fromstring(req.text)
        for isoform in uniprot:
            # get the sequence
            seq = isoform.find('{http://uniprot.org/uniprot}sequence')

            # get the accession number
            iso_accession = isoform.find('{http://uniprot.org/uniprot}accession')

            # add the values to the dictionary
            if seq.text and iso_accession.text:
                isoforms[iso_accession.text] = seq.text
        return isoforms
    except:
        AttributeError
        isoforms = {}
        return isoforms


def get_pdb_ids(protein_id):
    # Fetch PDB IDs associated with given UniProtID
    ssl._create_default_https_context = ssl._create_unverified_context
    url = 'https://www.uniprot.org/uploadlists/'
    params_ = {
        'from': 'ACC+ID',
        'to': 'PDB_ID',
        'format': 'tab',
        'query': protein_id.strip()
    }

    response = r.get(url, params=params_)
    response = response.text.split('\n')
    response = list(filter(None, response))

    pdbs = {}
    pdbs_per_protein = []
    for i in range(len(response)):
        try:
            if response[i].split('\t')[0] in pdbs.keys():
                pdbs_per_protein.append(response[i].split('\t')[1])
            elif response[i].split('\t')[0] not in pdbs.keys():
                pdbs_per_protein = []
                pdbs_per_protein.append(response[i].split('\t')[1])
        except:
            IndexError

        pdbs[response[i].split('\t')[0]] = pdbs_per_protein
    try:
        pdbs.pop('From')
    except:
        KeyError
        print(pdbs)
    return pdbs


def do_alignment(identifier, uniprotSequence, pdbSequence):
    print('Aligning ', identifier)
    if len(pdbSequence) >= 1:
        f = open(path_to_output_files + 'alignment_files/' + identifier + '_alignment' + ".txt",
                 "w")
        aligner.mode = 'local'
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        alignments = aligner.align(uniprotSequence, pdbSequence)
        alignments = (list(alignments))
        alignment_list = []
        for alignment in alignments:
            f.write(str(alignment))
            f.write('\n')
            f.write('\n')
            alignment = (str(alignment).strip().split('\n'))
            alignment = [''.join(['.' if m == ' ' else m for m in x]) for x in alignment]
            alignment_list.append(alignment)

    return alignment_list


def mutation_position_on_pdb(alignment_list, pos):
    which_alignment_to_go = 0
    for alignment in alignment_list:
        which_alignment_to_go += 1
        alignment_uniprot = alignment[0]
        alignment_pdb = alignment[2]
        startGap = 0
        if alignment_uniprot.startswith('.') or alignment_uniprot.startswith('-'):
            for k in alignment_uniprot:
                if k == '.' or k == '-':
                    startGap += 1
                else:
                    break
        countGap = startGap
        countResidue = 0
        canonicalRes = ' '
        pdbRes = ' '
        for j in alignment_uniprot[startGap:]:
            if j == '.' or j == '-':
                countGap += 1
            else:
                countResidue += 1
            if int(countResidue) == int(pos):
                canonicalRes = alignment_uniprot[countResidue + countGap - 1]
                try:  # Buraya sonradan try except ekledim.
                    pdbRes = alignment_pdb[countResidue + countGap - 1]
                except:
                    IndexError
                    pdbRes = 'nan'
                break
        if (alignment[1][countResidue + countGap - 1] == '|') or (alignment[1][countResidue + countGap - 1] == 'X'):
            if canonicalRes == pdbRes:
                pdb_alignStatus = 'aligned'
            elif canonicalRes != pdbRes:
                pdb_alignStatus = 'aligned*'
            countGap_pdb = 0
            countResidue_pdb = 0
            pdbRes = ' '
            for j in alignment_pdb[0:countResidue + countGap - 1]:
                if j == '.' or j == '-':
                    countGap_pdb += 1
            if alignment_pdb[countResidue + countGap - 1] == '.' or alignment_pdb[
                countResidue + countGap - 1] == '-':
                mutationPositionOnPDB = 'nan'  # there is no match, so position is Nan.
                posPDB = 'nan'
            else:
                posPDB = countResidue + countGap - countGap_pdb
                mutationPositionOnPDB = str(posPDB)
            break
        elif (canonicalRes == pdbRes) and ((alignment[1][countResidue + countGap - 1] == '.') or (
                alignment[1][countResidue + countGap - 1] == '-')):
            pdb_alignStatus = 'not_aligned'
            mutationPositionOnPDB = 'nan'
        elif (canonicalRes != pdbRes) and ((alignment[1][countResidue + countGap - 1] == '.') or (
                alignment[1][countResidue + countGap - 1] == '-')):
            pdb_alignStatus = 'not_aligned'
            mutationPositionOnPDB = 'nan'
        elif alignment_pdb[countResidue + countGap - 1] == '.' or alignment_pdb[
            countResidue + countGap - 1] == '-':
            mutationPositionOnPDB = 'nan'  # there is no match, so position is Nan.
            posPDB = 'nan'
    return (pdb_alignStatus, mutationPositionOnPDB, startGap, alignment_list[which_alignment_to_go - 1])


def find_position_on_pdb_for_range_annotations(posAnnotation, startGap, alignment_to_use):
    annotation_on_pdb_start = 'nan'
    annotation_on_pdb_end = 'nan'
    pos1 = int(posAnnotation.split('-')[0])
    count_gap = startGap
    count_residue = 0
    for j in alignment_to_use[0][startGap:]:
        if j == '.' or j == '-':
            count_gap += 1
        else:
            count_residue += 1
        if int(count_residue) == int(pos1):  # count gaps until the first position
            break
    annotation_on_up_start = int(pos1) + int(count_gap)

    pos2 = int(posAnnotation.split('-')[1])
    count_gap = startGap
    count_residue = 0
    for j in alignment_to_use[0][startGap:]:
        if j == '.' or j == '-':
            count_gap += 1
        else:
            count_residue += 1
        if int(count_residue) == int(pos2):  # count gaps until the first position
            break

    annotation_on_up_end = int(pos2) + int(count_gap)
    try:
        pdb_residue_start = alignment_to_use[2][annotation_on_up_start - 1].strip()
        if (pdb_residue_start == '.') or (pdb_residue_start == '-'):
            for ran in range(len(alignment_to_use[2][(annotation_on_up_start - 1):annotation_on_up_end])):
                if (alignment_to_use[2][(annotation_on_up_start - 1):annotation_on_up_end][ran] != '.') and \
                        (alignment_to_use[2][(annotation_on_up_start - 1):annotation_on_up_end][ran] != '-') and \
                        ((alignment_to_use[1][(annotation_on_up_start - 1):annotation_on_up_end][ran] == '|') or
                         (alignment_to_use[1][(annotation_on_up_start - 1):annotation_on_up_end][ran] == 'X')):
                    annotation_on_up_start += ran
                    break
        elif (pdb_residue_start != '.') and (pdb_residue_start != '-') and \
                ((alignment_to_use[1][annotation_on_up_start - 1] == '.') or (
                        alignment_to_use[1][annotation_on_up_start - 1] == '-')):
            for ran in range(len(alignment_to_use[2][(annotation_on_up_start - 1):annotation_on_up_end])):
                if ((alignment_to_use[1][(annotation_on_up_start - 1):annotation_on_up_end][ran] == '|') or
                        (alignment_to_use[1][(annotation_on_up_start - 1):annotation_on_up_end][ran] == 'X')):
                    annotation_on_up_start += ran
                    break
        count_gap_pdb = 0
        if annotation_on_up_start != 'nan':
            for q in alignment_to_use[2][0:annotation_on_up_start - 1]:
                if q == '.' or q == '-':
                    count_gap_pdb += 1
            if alignment_to_use[1][annotation_on_up_start] == '-' or alignment_to_use[1][annotation_on_up_start] == '.':
                annotation_on_pdb_start = 'nan'
            else:
                annotation_on_pdb_start = int(annotation_on_up_start) - count_gap_pdb
        else:
            annotation_on_pdb_start = 'nan'
    except:
        IndexError
    try:
        pdb_residue_end = alignment_to_use[2][annotation_on_up_end - 1].strip()
        if pdb_residue_end == '.' or pdb_residue_end == '-':
            for ran in range(len(alignment_to_use[2][(annotation_on_up_start - 1):annotation_on_up_end])):
                if ((alignment_to_use[1][annotation_on_up_start - 1:annotation_on_up_end][ran] == '.') or
                        (alignment_to_use[1][(annotation_on_up_start - 1):][ran] == '-')):
                    annotation_on_up_start += (ran - 1)
                    annotation_on_up_end = annotation_on_up_start
                    break
        elif (pdb_residue_end != '.') and (pdb_residue_end != '-') and \
                ((alignment_to_use[1][annotation_on_up_end - 1] == '.') or (
                        alignment_to_use[1][annotation_on_up_end - 1] == '-')):
            for ran in range(len(alignment_to_use[2][(annotation_on_up_start - 1):annotation_on_up_end])):
                if ((alignment_to_use[1][annotation_on_up_start - 1:annotation_on_up_end][ran] == '.') or
                        (alignment_to_use[1][(annotation_on_up_start - 1):][ran] == '-')):
                    annotation_on_up_start += (ran - 1)
                    annotation_on_up_end = annotation_on_up_start
                    break
        count_gap_pdb = 0
        if annotation_on_up_end != 'nan':
            for q in alignment_to_use[2][0:annotation_on_up_end - 1]:
                if q == '.' or q == '-':
                    count_gap_pdb += 1
            if alignment_to_use[1][annotation_on_up_end - 1] == '-' or alignment_to_use[1][
                annotation_on_up_end - 1] == '.' and annotation_on_pdb_start == 'nan':
                annotation_on_pdb_end = 'nan'
            elif alignment_to_use[1][annotation_on_up_end - 1] == '-' or alignment_to_use[1][
                annotation_on_up_end - 1] == '.' and annotation_on_pdb_start != 'nan':
                annotation_on_pdb_end = int(annotation_on_up_end) - count_gap_pdb
            else:
                annotation_on_pdb_end = int(annotation_on_up_end) - count_gap_pdb
        else:
            annotation_on_pdb_end = 'nan'
    except:
        IndexError  # Say isoform 2 is matched with the length 100, but canonical is 150 aa long. If there is an annotation at 105. position, for the isoform it throws an index error.

    if annotation_on_pdb_start == 'nan' and annotation_on_pdb_end != 'nan':
        annotation_on_pdb_start = annotation_on_up_start - count_gap_pdb
    if annotation_on_pdb_start == annotation_on_pdb_end:
        annotation_on_pdb_start = 'nan'
        annotation_on_pdb_end = 'nan'
    return annotation_on_up_start, annotation_on_up_end, annotation_on_pdb_start, annotation_on_pdb_end


def annotation_pos_on_pdb(annot_positions, startGap, alignment_to_use, identifier):
    newpos = []
    if annot_positions != 'nan':  # if not Nan #We need to find up_gap for each annotation position.
        annot_positions = (str(annot_positions).replace("'", ''))
        annot_positions = (str(annot_positions).replace('[', ''))
        annot_positions = (str(annot_positions).replace("]", ''))
        positionList_perAnnotation = annot_positions.split(',')
        positionList_perAnnotation = [h.strip() for h in positionList_perAnnotation]
        position_start_on_pdb = 'nan'
        position_end_on_pdb = 'nan'

        for position in range(len(positionList_perAnnotation)):
            if ('-' not in str(positionList_perAnnotation[position])) and (
                    str(positionList_perAnnotation[position]) != '?') and (
                    str(positionList_perAnnotation[position]) != '') and (
                    len(str(positionList_perAnnotation[position])) != 0):
                count_gap = startGap
                count_residue = 0
                for j in alignment_to_use[0][startGap:]:
                    if j == '.' or j == '-':
                        count_gap += 1
                    else:
                        count_residue += 1
                    try:
                        if int(count_residue) == int(positionList_perAnnotation[position]):
                            break
                    except:
                        ValueError

                annotation_on_up = int(positionList_perAnnotation[position]) + int(count_gap)
                try:
                    pdb_residue_start = alignment_to_use[2][annotation_on_up - 1].strip()
                except:
                    IndexError
                    pdb_residue_start = 'nan'
                if pdb_residue_start != 'nan':
                    try:
                        if (pdb_residue_start == '.') or (pdb_residue_start == '-'):
                            for ran in range(len(alignment_to_use[2][(annotation_on_up - 1):annotation_on_up])):
                                if (alignment_to_use[2][(annotation_on_up - 1):annotation_on_up][ran] != '.') and \
                                        (alignment_to_use[2][(annotation_on_up - 1):annotation_on_up][
                                             ran] != '-') and \
                                        ((alignment_to_use[1][(annotation_on_up - 1):annotation_on_up][
                                              ran] == '|') or
                                         (alignment_to_use[1][(annotation_on_up - 1):annotation_on_up][
                                              ran] == 'X')):
                                    annotation_on_up += ran
                                    break
                        elif (pdb_residue_start != '.') and (pdb_residue_start != '-') and \
                                ((alignment_to_use[1][annotation_on_up - 1] == '.') or (
                                        alignment_to_use[1][annotation_on_up - 1] == '-')):
                            for ran in range(len(alignment_to_use[2][(annotation_on_up - 1):annotation_on_up])):
                                if ((alignment_to_use[1][(annotation_on_up - 1):annotation_on_up][ran] == '|') or
                                        (alignment_to_use[1][(annotation_on_up - 1):annotation_on_up][ran] == 'X')):
                                    annotation_on_up += ran
                                    break
                        count_gap_pdb = 0
                        for q in alignment_to_use[2][0:annotation_on_up - 1]:
                            if q == '.' or q == '-':
                                count_gap_pdb += 1
                        if alignment_to_use[1][annotation_on_up] == '-' or alignment_to_use[1][
                            annotation_on_up] == '.':
                            annotation_on_pdb = 'nan'
                        else:
                            annotation_on_pdb = int(annotation_on_up) - count_gap_pdb
                        # annotation_on_pdb = int(annotation_on_up) - count_gap_pdb
                        if count_gap_pdb == annotation_on_up:
                            annotation_on_pdb = 'nan'
                        try:
                            if alignment_to_use[2][count_gap_pdb + annotation_on_pdb - 1] == '.' or alignment_to_use[2][
                                count_gap_pdb + annotation_on_pdb - 1] == '-':
                                annotation_on_pdb = 'nan'
                        except:
                            IndexError
                            annotation_on_pdb = 'nan'
                    except:
                        IndexError
                        annotation_on_pdb = 'nan'

                    newpos.append(annotation_on_pdb)

            elif ('-' in str(positionList_perAnnotation[position])) and (
                    str(positionList_perAnnotation[position]) != '?') and (
                    str(positionList_perAnnotation[position]) != ' ') and (
                    len(str(positionList_perAnnotation[position])) != 0):
                try:
                    position_start_on_pdb = \
                        find_position_on_pdb_for_range_annotations(positionList_perAnnotation[position],
                                                                   startGap, alignment_to_use)[2]
                    position_end_on_pdb = \
                        find_position_on_pdb_for_range_annotations(positionList_perAnnotation[position],
                                                                   startGap, alignment_to_use)[3]
                except:
                    ValueError
                newpositions = str(position_start_on_pdb) + '-' + str(position_end_on_pdb)
                newpos.append(newpositions)  # positions on the pdb Fasta sequence.
            else:
                pass
                # print('positionList_perAnnotation[position]', positionList_perAnnotation[position], identifier)
    return newpos


def final_stage(df):
    for i in df.index:
        print(i, 'of', len(df))
        identifier = df.at[i, 'uniprotID'] + '_' + df.at[i, 'pdbID'] + '_' + df.at[i, 'chain'] + '_'
        alignment_list = do_alignment(identifier, df.at[i, 'uniprotSequence'], df.at[i, 'pdbSequence'])
        df.at[i, 'pdb_alignStatus'] = mutation_position_on_pdb(alignment_list, df.at[i, 'pos'])[0]
        df.at[i, 'mutationPositionOnPDB'] = mutation_position_on_pdb(alignment_list, df.at[i, 'pos'])[1]
        startGap = mutation_position_on_pdb(alignment_list, df.at[i, 'pos'])[2]
        alignment_to_use = mutation_position_on_pdb(alignment_list, df.at[i, 'pos'])[3]
        for annot in annotation_list:
            df.at[i, annot] = annotation_pos_on_pdb(df.at[i, annot], startGap, alignment_to_use, identifier)
        if str(df.at[i, 'domStart']) != 'nan' and str(df.at[i, 'domEnd']) != 'nan' and \
                ((str(df.at[i, 'domStart']) != '-1' and str(df.at[i, 'domEnd']) != '-1' and
                  str(df.at[i, 'domStart']) != '-1.0' and str(df.at[i, 'domEnd']) != '-1.0')):
            domainLoc = str(df.at[i, 'domStart']).split('.')[0] + '-' + str(df.at[i, 'domEnd']).split('.')[0]
            domain_pos = find_position_on_pdb_for_range_annotations(domainLoc, startGap, alignment_to_use)
            df.at[i, 'domainStartonPDB'] = domain_pos[2]
            df.at[i, 'domainEndonPDB'] = domain_pos[3]
        elif str(df.at[i, 'domStart']) != '-1' or str(df.at[i, 'domEnd']) != '-1' or \
                str(df.at[i, 'domStart']) != '-1.0' or str(df.at[i, 'domEnd']) != '-1.0':
            df.at[i, 'domainStartonPDB'] = 'nan'
            df.at[i, 'domainEndonPDB'] = 'nan'
    return df


def alignment(dataframe_to_align):
    domainList = ['domStart', 'domEnd']
    result = final_stage(dataframe_to_align)
    return result


def run_freesasa(infile, outfile, include_hetatms=True, outdir=None, force_rerun=False):
    if not outdir:
        outdir = ''

    outfile = op.join(outdir, outfile)
    if ssbio.utils.force_rerun(flag=force_rerun, outfile=outfile):
        if include_hetatms:
            shell_command = 'freesasa --format=rsa --hetatm {} -o {}'.format(infile, outfile)
        else:
            shell_command = 'freesasa --format=rsa {} -o {}'.format(infile, outfile)
        command = subprocess.Popen(shell_command,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   shell=True)
        out, err = command.communicate()

    return outfile


def distance(x1, y1, z1, x2, y2, z2):
    d = math.sqrt(math.pow(x2 - x1, 2) +
                  math.pow(y2 - y1, 2) +
                  math.pow(z2 - z1, 2) * 1.0)
    return d


def find_distance(coordMut, coordAnnot):
    try:
        dist = distance(float(coordMut[0]), float(coordMut[1]), float(coordMut[2]), float(coordAnnot[0]),
                        float(coordAnnot[1]), float(coordAnnot[2]))
        return "%.2f" % dist
    except:  # mesela pdbdeki residue seqreste görünemedi kısaydı, koordinat yok o yüzden.
        ValueError
        dist = 'nan'
        return dist


def threeToOne(variant):
    if variant == "ALA":
        variant = "A"
    elif variant == "ARG":
        variant = "R"
    elif variant == "VAL":
        variant = "V"
    elif variant == "GLU":
        variant = "E"
    elif variant == "PRO":
        variant = "P"
    elif variant == "LEU":
        variant = "L"
    elif variant == "GLY":
        variant = "G"
    elif variant == "ASN":
        variant = "N"
    elif variant == "SER":
        variant = "S"
    elif variant == "GLN":
        variant = "Q"
    elif variant == "THR":
        variant = "T"
    elif variant == "MET":
        variant = "M"
    elif variant == "LYS":
        variant = "K"
    elif variant == "ASP":
        variant = "D"
    elif variant == "ILE":
        variant = "I"
    elif variant == "PHE":
        variant = "F"
    elif variant == "TRP":
        variant = "W"
    elif variant == "TYR":
        variant = "Y"
    elif variant == "HIS":
        variant = "H"
    elif variant == "CYS":
        variant = "C"
    elif variant == 'UNK':
        variant = 'X'
    elif variant == 'ASX':
        variant = 'O'
    return (variant)


def get_coords(annot, alignments):
    for alignment in alignments[0]:
        alignment = (str(alignment).strip().split('\n'))
        startGap = 0
        if alignment[0].startswith('.'):
            for k in alignment[0]:
                if k == '.' or k == '-':
                    startGap += 1
                else:
                    break
        countGap = startGap
        countResidue = 0
        for j in alignment[0][startGap:]:
            if j == '.' or j == '-':
                countGap += 1
            else:
                countResidue += 1
            if countResidue == float(annot):  # mutPDB pos on pdb mesela, ya da diğer annor positions
                break
        countGap_pdb = 0
        countResidue_pdb = 0
        for m in alignment[2][0:countResidue + countGap - 1]:
            if m == '.' or m == '-':
                countGap_pdb += 1
        posAtom = countResidue + countGap - countGap_pdb

        realpdbStart = 0
        for j in alignment[2]:
            if j == '.' or j == '-':
                realpdbStart += 1
            else:
                break

        if (alignment[2][countResidue + countGap - 1] != '-') and (float(annot) >= float(realpdbStart) + 1):
            try:
                coordinates = alignments[1]
                residue_numbers = alignments[2]
                coordWeWant = coordinates[posAtom - 1]
                residue_number_we_want = residue_numbers[posAtom - 1]

            except:
                IndexError
                coordWeWant = 'nan'
        else:
            coordWeWant = 'nan'
        return coordWeWant, posAtom, residue_number_we_want


def get_alignments_3D(data, pdb_path):  # input sayı; output coordinates
    for name in glob.glob(pdb_path):
        atomSequence = ''
        coords = []
        resnums_for_sasa = []
        with open(name, encoding="utf8") as f:
            for line in f.readlines():
                if data.at[i, 'source'] != 'MODBASE':
                    if line[0:4].strip() == 'ATOM' and line[13:15].strip() == 'CA' and line[21].upper() == data.at[
                        i, 'chain'].upper():
                        atomSequence += threeToOne(line[17:20].strip())
                        coords.append([line[31:38].strip(), line[39:46].strip(), line[47:54].strip()])
                        resnums_for_sasa.append(line[22:26].strip())
                    elif line[0:4].strip() == 'ATOM' and line[13:15].strip() == 'CA' and line[21] == ' ':
                        atomSequence += threeToOne(line[17:20].strip())
                        coords.append([line[31:38].strip(), line[39:46].strip(), line[47:54].strip()])
                        resnums_for_sasa.append(line[22:26].strip())
                else:
                    if line[0:7].strip() == 'ATOM' and line[13:15].strip() == 'CA':
                        atomSequence += threeToOne(line[17:20].strip())
                        coords.append([line[31:38].strip(), line[39:46].strip(), line[47:54].strip()])
                        resnums_for_sasa.append(line[22:26].strip())

        f = open(path_to_output_files + '3D_alignment/' + data.at[
            i, 'uniprotID'] + '_' + data.at[i, 'pdbID'] + '_' + str(data.at[i, 'chain']) + '_alignment' + ".txt",
                 "w")

        aligner.mode = 'local'
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        # print('pdb', data.at[i, 'pdbSequence'],len(data.at[i, 'pdbSequence']), type(data.at[i, 'pdbSequence']))
        # print('atomSeq', atomSequence, len(atomSequence), type(atomSequence))
        # print(data.at[
        #    i, 'uniprotID'] + '_' + data.at[i, 'pdbID'] + '_' + str(data.at[i, 'chain']) )
        alignments = aligner.align(data.at[i, 'pdbSequence'], atomSequence)
        alignments = (list(alignments))
        for alignment in alignments:
            f.write(str(alignment))
            f.write('\n')
            f.write('\n')
        return alignments, coords, resnums_for_sasa


def sasa(data, sasa_pos):
    sasa = 'nan'
    for filename in glob.glob(path_to_output_files + 'freesasa_files/*'):
        if data.at[i, 'source'] == 'PDB':
            fname = filename.split('.')[0].split('/')[-1].upper()
        elif data.at[i, 'source'] == 'MODBASE':
            fname = filename.split('.')[0].split('/')[-1]
        elif data.at[i, 'source'] == 'SWISSSMODEL':
            fname = filename.split('_')[2]

        if data.at[i, 'pdbID'] == fname:
            files = open(filename, 'r')
            file = files.readlines()
            for k in file:
                if k.strip()[10:13] == sasa_pos:
                    residue = str(k[4:7].strip())
                    if data.at[i, 'wt'] == threeToOne(residue):
                        sasa = str(k[22:28]).strip('\n')
                        return (sasa)
                    elif data.at[i, 'wt'] != threeToOne(residue):
                        sasa = str(k[22:28]).strip('\n') + '*'
                        return (sasa)
                    else:
                        return 'nan'  #######


def standardize(df, get_columns):
    cols_to_change = ['sasa', 'domaindistance3D', 'disulfide', 'intMet', 'intramembrane',
                      'naturalVariant', 'dnaBinding', 'activeSite', 'nucleotideBinding',
                      'lipidation', 'site', 'transmembrane', 'crosslink', 'mutagenesis',
                      'strand', 'helix', 'turn', 'metalBinding', 'repeat', 'caBinding',
                      'topologicalDomain', 'bindingSite', 'region', 'signalPeptide',
                      'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil', 'peptide',
                      'transitPeptide', 'glycosylation', 'propeptide']
    for col in cols_to_change:  # because in the other ones, they are 3D distance. Here, no distance calculated.
        df[col] = 'nan'
    df = df[get_columns.columns]

    return df


def get_interface_positions(dataframe, column1, column2):
    interface_positions = {}
    for i in dataframe.index:
        if dataframe.at[i, column1] not in interface_positions and dataframe.at[i, column1 + '_IRES'] != '[]':
            interface_positions[dataframe.at[i, column1]] = dataframe.at[i, str(column1 + '_IRES')]
        elif dataframe.at[i, column1] in interface_positions and dataframe.at[i, column1 + '_IRES'] != '[]':
            interface_positions[dataframe.at[i, column1]] = interface_positions[dataframe.at[i, column1]].strip(
                ']') + ',' + (dataframe.at[i, str(column1 + '_IRES')]).strip('[')
        if dataframe.at[i, column2] not in interface_positions and dataframe.at[i, column2 + '_IRES'] != '[]':
            interface_positions[dataframe.at[i, column2]] = dataframe.at[i, str(column2 + '_IRES')]
        elif dataframe.at[i, column2] in interface_positions and dataframe.at[i, column2 + '_IRES'] != '[]':
            interface_positions[dataframe.at[i, column2]] = interface_positions[dataframe.at[i, column2]].strip(
                ']') + ',' + (dataframe.at[i, str(column2 + '_IRES')]).strip('[')

    try:
        for key, value in interface_positions.items():
            n = []
            m = []
            if value != '[]':
                valueList = value.split(',')
                valueList[0] = str(valueList[0]).strip('[')
                valueList[-1] = str(valueList[-1]).strip(']')
                for val in valueList:
                    if '-' in val:
                        for r in range(int(val.split('-')[0]), int(val.split('-')[1]) + 1):
                            n.append(r)
                    else:
                        m.append(int(val))
                fin = m + n

                interface_positions[key] = fin
    except:
        ValueError

    return interface_positions


aligner = Align.PairwiseAligner()
"""
STEP 1
Get input data as a console input.
Add datapoint identifier and remove non-standard input.
"""
input_set = input('Enter Query DataPoint\n'
                  'Option 1: Comma-separated list of idenfiers (UniProt ID-wt residue-position-mutated residue (e.g. Q9Y4W6-N-432-T or Q9Y4W6-N-432-T, Q9Y4W6-N-432-T)) \n'
                  'Option 2: Enter comma-separated file path ')

path_to_input_files = 'input_files/'
path_to_output_files = 'out_files/'
path_to_domains = path_to_input_files + 'domains.txt'
swiss_model_path = path_to_input_files + 'INDEX.json'
fisher_path = path_to_input_files + 'significant_domains.txt'
path_to_interfaces = path_to_input_files + 'H_sapiens_interfacesHQ.txt'
outpath = path_to_output_files + 'feature_vector.txt'
buffer = path_to_output_files + 'file_so_far.txt'


sys.stdout = open('log.txt', 'w')
print('Creating directories...')
os.makedirs(path_to_output_files, exist_ok=True)
os.makedirs(path_to_output_files + 'pdb_structures/', exist_ok=True)
os.makedirs(path_to_output_files + 'alignment_files/', exist_ok=True)
os.makedirs(path_to_output_files + 'swissmodel_structures/', exist_ok=True)
os.makedirs(path_to_output_files + 'modbase_structures/', exist_ok=True)
os.makedirs(path_to_output_files + 'modbase_structures_individual/', exist_ok=True)
os.makedirs(path_to_output_files + 'freesasa_files/', exist_ok=True)
os.makedirs(path_to_output_files + '3D_alignment/', exist_ok=True)


annotation_list = ['disulfide', 'intMet', 'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                   'nucleotideBinding', 'lipidation', 'site', 'transmembrane', 'crosslink', 'mutagenesis', 'strand',
                   'helix', 'turn', 'metalBinding', 'repeat', 'topologicalDomain', 'caBinding', 'bindingSite', 'region',
                   'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil', 'peptide',
                   'transitPeptide', 'glycosylation', 'propeptide']

print('Feature vector generation started...\n')

data = pd.DataFrame()

try:
    if '-' in input_set:
        data = data.append(pd.Series([j.strip() for j in input_set.split('-')]), ignore_index=True)
        data.columns = ['uniprotID', 'wt', 'pos', 'mut']
    if ',' in input_set:
        input_set = [i.strip() for i in input_set.split(',')]
        for i in input_set:
            data = data.append(pd.Series([j.strip() for j in i.split('-')]), ignore_index=True)
        data.columns = ['uniprotID', 'wt', 'pos', 'mut']
    elif '\t' in input_set:
        input_set = [i.strip() for i in input_set.split('\t')]
        for i in input_set:
            data = data.append(pd.Series([j.strip() for j in i.split('-')]), ignore_index=True)
        data.columns = ['uniprotID', 'wt', 'pos', 'mut']
    elif ' ' in input_set:
        input_set = [i.strip() for i in input_set.split(' ')]
        for i in input_set:
            data = data.append(pd.Series([j.strip() for j in i.split('-')]), ignore_index=True)
        data.columns = ['uniprotID', 'wt', 'pos', 'mut']
    elif '.txt' in input_set:
        data = pd.read_csv(input_set, sep='\t', names=['uniprotID', 'wt', 'pos', 'mut'])
except:
    ValueError
    print('Please check the input format.')

# Exclude termination codons, synonymous mutations and any non-standard residues such as Sec, 4 or 6.
aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
data.wt = data.wt.str.strip()
data.mut = data.mut.str.strip()
data = data[data.wt.isin(aa_list)]
data = data[data.mut.isin(aa_list)]

for i in data.index:
    data.at[i, 'datapoint'] = data.at[i, 'uniprotID'] + data.at[i, 'wt'] + str(data.at[i, 'pos']) + data.at[i, 'mut']

data = data.astype(str)
data = data[data.uniprotID != 'Q8WZ42']

if len(data) == 0:
    print('Feature vectore generation terminated.')
else:
    """
    STEP 2
    Add physicochemical properties.
    """
    print('Adding physicochemical properties...\n')

    data = add_physicochemical(data)

    """
    STEP 3
    Add domains-related information.
    """
    print('Adding domains\n')

    data = add_domains(data, path_to_domains)

    data = data.astype(str)
    data = data.replace({'NaN': 'nan'})
    data.domain = data.domain.replace({'nan': '-1'})  # Domainler nan ise -1 olsun.
    data.domStart = data.domStart.replace({'nan': '-1'})
    data.domEnd = data.domEnd.replace({'nan': '-1'})
    data.distance = data.distance.replace({'nan': '-1'})

    """
    STEP 4
    Retrieve canonical and isoform UniProt sequences.
    Add to the data frame.
    """
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
    data = data.astype(str)
    data['whichIsoform'] = 'nan'
    data.replace({'': 'nan'}, inplace=True)
    data['wt_sequence_match'] = ''
    for i in data.index:
        if len(data.at[i, 'uniprotSequence']) >= int(data.at[i, 'pos']):
            wt = data.at[i, 'wt']
            can = str(data.at[i, 'uniprotSequence'])[int(data.at[i, 'pos']) - 1]
            if wt == can:
                data.at[i, 'wt_sequence_match'] = 'm'
            elif wt != can:
                isoList = isoform_fasta[isoform_fasta['uniprotID'] == data.at[i, 'uniprotID']].isoformSequence.to_list()
                for k in isoList:
                    if len(k) >= int(data.at[i, 'pos']):
                        resInIso = k[int(int(data.at[i, 'pos']) - 1)]
                        if wt == resInIso:
                            whichIsoform = isoform_fasta[isoform_fasta.isoformSequence == k].whichIsoform.to_list()[0]
                            data.at[i, 'wt_sequence_match'] = 'i'
                            data.at[i, 'whichIsoform'] = whichIsoform
                            break

        elif len(data.at[i, 'uniprotSequence']) < int(data.at[i, 'pos']):
            isoList = isoform_fasta[isoform_fasta['uniprotID'] == data.at[i, 'uniprotID']].isoformSequence.to_list()
            for k in isoList:
                if len(k) >= int(data.at[i, 'pos']):
                    resInIso = k[int(int(data.at[i, 'pos']) - 1)]
                    wt = data.at[i, 'wt']
                    if wt == resInIso:
                        whichIsoform = isoform_fasta[isoform_fasta.isoformSequence == k].whichIsoform.to_list()[0]
                        data.at[i, 'wt_sequence_match'] = 'i'
                        data.at[i, 'whichIsoform'] = whichIsoform
                        break

    data.wt_sequence_match = data.wt_sequence_match.astype('str')
    data.replace({'': 'nan'}, inplace=True)
    data_size = len(data.drop_duplicates(['datapoint']))
    not_match_in_uniprot = data[(data.uniprotSequence == 'nan') | (data.wt_sequence_match == 'nan')]
    uniprot_matched = data[(data.uniprotSequence != 'nan') & (data.wt_sequence_match != 'nan')]
    data = None

    print('You have %d data points that failed to match a UniProt Sequence\nProceeding with %d remaining...\n'
          % (len(not_match_in_uniprot.drop_duplicates(['datapoint'])),
             len(uniprot_matched.drop_duplicates(['datapoint']))))

    """
    STEP 5
    Retrieve related PDB sequences, extract their sequences.
    Add to the data frame.
    """

    pdb_fasta = pd.DataFrame(columns=['pdbID', 'chain', 'pdbSequence'])
    pdb_info = pd.DataFrame(columns=['uniprotID', 'pdbID', 'chain', 'resolution'])

    print('Retrieving PDB structures...\n')
    pdbs = []
    for protein in list(set(uniprot_matched.uniprotID.to_list())):
        dict = get_pdb_ids(protein)
        pdbs.append([dict[k] for k in dict])

    print('Processing PDB structures...\n')
    try:
        pdbs = [j.strip('[').strip(']').strip().strip('\'').strip('\"') for j in
                ((',').join([str(item) for item in pdbs])).split(',')]
    except IndexError:
        pdbs = []
        print('No PDB structure found for the query. ')

    pdbs = list(filter(None, pdbs))
    pdbs = (set(pdbs))
    pdbs = [i.lower() for i in pdbs]
    pdbl = PDBList()
    parser = PDBParser()
    index = 0
    print('Starting PDB structures download...\n')
    import shutil
    try:
        shutil.rmtree('obsolete')
    except OSError as e:
        pass
    existing_pdb = glob.glob(path_to_output_files + 'pdb_structures/*')
    existing_pdb = [i.split('/')[-1].split('.')[0].lower() for i in existing_pdb]
    cnt = 0
    for search in pdbs:
        print(cnt, 'of', len(pdbs))
        try:
            if search.lower() not in existing_pdb:
                file = pdbl.retrieve_pdb_file(search, pdir=path_to_output_files + 'pdb_structures/', file_format="pdb")
            else:
                print('PDB structure file exists..')
                for filename in os.listdir(path_to_output_files + 'pdb_structures/'):
                    os.rename(path_to_output_files + 'pdb_structures/' + filename,
                                  path_to_output_files + 'pdb_structures/' + ''.join(filename.split('.')[0])+ '.pdb')

                file = path_to_output_files + 'pdb_structures/' + search + '.pdb'

                base = os.path.splitext(file)[0]
                base = '/'.join(base.split('/')[0:-1]) + '/pdb' + base.split('/')[-1]
                os.rename(file, base + ".ent")
                file = base + '.ent'

            resolution_method = parser.get_structure(search, file)
            for record in SeqIO.parse(file, "pdb-seqres"):
                if record.dbxrefs[0].split(':')[0] == 'UNP':
                    pdb_fasta.at[index, 'pdbID'] = record.id.split(':')[0]
                    pdb_fasta.at[index, 'chain'] = record.id.split(':')[1]
                    pdb_fasta.at[index, 'pdbSequence'] = str(record.seq)
                    pdb_info.at[index, 'uniprotID'] = record.dbxrefs[0].split(':')[1]
                    pdb_info.at[index, 'pdbID'] = record.id.split(':')[0]
                    pdb_info.at[index, 'chain'] = record.annotations["chain"]
                    pdb_info.at[index, 'resolution'] = resolution_method.header['resolution']
                index += 1
        except:
            IndexError
            pdb_info.at[index, 'uniprotID'] = 'nan'
            pdb_info.at[index, 'pdbID'] = 'nan'
            pdb_info.at[index, 'chain'] = 'nan'
            pdb_info.at[index, 'resolution'] = 'nan'
        cnt +=1
    print()
    print('PDB file processing finished..')
    for filename in os.listdir(path_to_output_files + 'pdb_structures/'):
        try:
            os.rename(path_to_output_files + 'pdb_structures/' + filename,
                      path_to_output_files + 'pdb_structures/' + filename[:-4] + '.txt')
        except:
            FileNotFoundError

    for filename in os.listdir(path_to_output_files + 'pdb_structures/'):
        try:
            if filename.startswith("pdb"):
                os.rename(path_to_output_files + 'pdb_structures/' + filename,
                          path_to_output_files + 'pdb_structures/' + filename[3:])
        except:
            FileNotFoundError
    # pdb_info.to_csv('/Users/fatmacankara/Desktop/new_benchmark/pdb_info.txt', sep='\t', index=False)
    # pdb_fasta.to_csv('/Users/fatmacankara/Desktop/new_benchmark/pdb_fasta.txt', sep='\t', index=False)

    uniprot_matched = pd.merge(uniprot_matched, pdb_info, on='uniprotID', how='left')
    uniprot_matched = uniprot_matched.astype(str)
    uniprot_matched = uniprot_matched.drop_duplicates()

    uniprot_matched = uniprot_matched.merge(pdb_fasta, on=['pdbID', 'chain'], how='left')
    uniprot_matched = uniprot_matched.astype(str)

    with_pdb = uniprot_matched[(uniprot_matched.pdbID != 'nan') & (
            (uniprot_matched.resolution != 'nan') & (uniprot_matched.resolution != 'OT') & (
            uniprot_matched.resolution != 'None'))].drop_duplicates()
    no_pdb = uniprot_matched[(uniprot_matched.pdbID == 'nan') | (
            (uniprot_matched.resolution == 'nan') | (uniprot_matched.resolution == 'OT') | (
            uniprot_matched.resolution == 'None'))]
    no_pdb = no_pdb[~no_pdb.datapoint.isin(with_pdb.datapoint.to_list())]
    no_pdb.drop(columns=['chain', 'pdbID', 'pdbSequence', 'resolution'], inplace=True)

    print(
        'PDB Information successfully added...\nPDB structures are found for %d of %d.\n%d of %d failed to match with PDB structure.\n'
        % (len(with_pdb.drop_duplicates(['datapoint'])), len(uniprot_matched.drop_duplicates(['datapoint'])),
           len(no_pdb.drop_duplicates(['datapoint'])), len(uniprot_matched.drop_duplicates(['datapoint']))))

    with_pdb = with_pdb.sort_values(['uniprotID', 'resolution'], axis=0, ascending=True)
    with_pdb = with_pdb.drop_duplicates(['uniprotID', 'wt', 'mut', 'pos', 'pdbSequence'], keep='first')
    with_pdb.replace({'': 'nan'}, inplace=True)

    if len(with_pdb) == 0:
        with_pdb['pdbInfo'] = ''
    else:
        for i in with_pdb.index:
            try:
                res = str(with_pdb.at[i, 'resolution'])
                chain = with_pdb.at[i, 'chain']
                new = with_pdb.at[i, 'pdbID'] + ':' + chain + ':' + res
                with_pdb.at[i, 'pdbInfo'] = new
            except:
                TypeError
                with_pdb.at[i, 'pdbInfo'] = 'nan'

    with_pdb = with_pdb[['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                         'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence', 'pdbSequence',
                         'wt_sequence_match',
                         'whichIsoform', 'pdbID', 'resolution', 'chain', 'pdbInfo', 'datapoint']]

    # with_pdb.to_csv('/Users/fatmacankara/Desktop/new_benchmark/with_pdb.txt', sep='\t', index=False)
    # no_pdb.to_csv('/Users/fatmacankara/Desktop/new_benchmark/no_pdb.txt', sep='\t', index=False)

    # If the query data points are found in no_match_in_uniprot data frame, it will not give any results.
    # If the query data points are found in no_pdb data frame, it will be searched in the modbase and swiss_model steps.
    # If the query data points are found in with_pdb data frame, it will be searched in the following steps.

    """
    STEP 6
    Retrieve sequence annotations.
    Add to the data frame.
    """

    if len(with_pdb) > 0:
        with_pdb = add_annotations(with_pdb)
    else:
        new_cols = with_pdb.columns.to_list() + ['disulfide', 'intMet', 'intramembrane', 'naturalVariant', 'dnaBinding',
                                                 'activeSite',
                                                 'nucleotideBinding', 'lipidation', 'site', 'transmembrane',
                                                 'crosslink', 'mutagenesis', 'strand',
                                                 'helix', 'turn', 'metalBinding', 'repeat', 'topologicalDomain',
                                                 'caBinding', 'bindingSite', 'region',
                                                 'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif',
                                                 'coiledCoil', 'peptide',
                                                 'transitPeptide', 'glycosylation', 'propeptide', 'disulfideBinary',
                                                 'intMetBinary', 'intramembraneBinary',
                                                 'naturalVariantBinary', 'dnaBindingBinary', 'activeSiteBinary',
                                                 'nucleotideBindingBinary', 'lipidationBinary', 'siteBinary',
                                                 'transmembraneBinary', 'crosslinkBinary', 'mutagenesisBinary',
                                                 'strandBinary', 'helixBinary', 'turnBinary', 'metalBindingBinary',
                                                 'repeatBinary', 'topologicalDomainBinary', 'caBindingBinary',
                                                 'bindingSiteBinary', 'regionBinary', 'signalPeptideBinary',
                                                 'modifiedResidueBinary', 'zincFingerBinary', 'motifBinary',
                                                 'coiledCoilBinary', 'peptideBinary', 'transitPeptideBinary',
                                                 'glycosylationBinary', 'propeptideBinary']
        with_pdb = pd.DataFrame(columns=new_cols)
    try:
        with_pdb.whichIsoform = with_pdb.whichIsoform.astype('str')
    except:
        AttributeError
        with_pdb['whichIsoform'] = ''

    with_pdb = with_pdb.astype(str)
    with_pdb = with_pdb.replace({'NaN': 'nan'})
    with_pdb.replace({'[]': 'nan'}, inplace=True)
    with_pdb.replace({'nan-nan': 'nan'}, inplace=True)
    with_pdb.replace({'': 'nan'}, inplace=True)

    # with_pdb.to_csv('/Users/fatmacankara/Desktop/new_benchmark/with_pdb_with_annotations.txt', sep='\t', index=False)
    """
    STEP 7
    Do alignment for PDB
    """
    # Canonical matches, i.e. labelled as m, canonical sequences will be aligned with PDB sequences.
    # Isoform matches, i.e. labelled as i, isoform sequences will be aligned with PDB sequences.
    with_pdb['uniprotSequence'] = with_pdb['uniprotSequence'].str.replace('U', 'C')
    with_pdb['pdbSequence'] = with_pdb['pdbSequence'].str.replace('U', 'C')

    dfM = with_pdb[with_pdb.wt_sequence_match == 'm']
    dfM = dfM.sort_values(['uniprotID', 'resolution'], axis=0, ascending=True)
    dfM = dfM.drop_duplicates(['uniprotID', 'wt', 'mut', 'pos', 'pdbSequence'], keep='first')

    dfNM = with_pdb[with_pdb.wt_sequence_match == 'i']
    dfNM = dfNM.sort_values(['uniprotID', 'resolution'], axis=0, ascending=True)
    dfNM = dfNM.drop_duplicates(['uniprotID', 'wt', 'mut', 'pos', 'pdbSequence'], keep='first')
    dfNM.rename(columns={'isoformSequence': 'uniprotSequence'}, inplace=True)

    dfM = dfM.astype(str)
    dfNM = dfNM.astype(str)

    dfM.reset_index(inplace=True)
    dfM.drop(['index'], axis=1, inplace=True)
    dfNM.reset_index(inplace=True)
    dfNM.drop(['index'], axis=1, inplace=True)
    # dfM.to_csv('/Users/fatmacankara/Desktop/new_benchmark/dfM.txt', sep='\t', index=False)
    # dfNM.to_csv('/Users/fatmacankara/Desktop/new_benchmark/dfNM.txt', sep='\t', index=False)
    uniprot_matched_size = len(uniprot_matched.drop_duplicates(['datapoint']))
    uniprot_matched = None
    pdb_fasta = None
    pdb_info = None
    pdbs = None
    existing_pdb = None
    with_pdb_size = len(with_pdb.drop_duplicates(['datapoint']))
    with_pdb = None


    print('Aligning sequences...\n')

    aligned_m = final_stage(dfM)
    aligned_nm = final_stage(dfNM)

    # When PDB sequence is nan, it is wrongly aligned to the UniProt sequence. Fix them.
    for i in aligned_m.index:
        if aligned_m.at[i, 'pdbSequence'] == 'nan':
            aligned_m.at[i, 'mutationPositionOnPDB'] = 'nan'
            aligned_m.at[i, 'domainStartonPDB'] = 'nan'
            aligned_m.at[i, 'domainEndonPDB'] = 'nan'
            aligned_m.at[i, 'pdb_alignStatus'] = 'nan'

    for i in aligned_nm.index:
        if aligned_nm.at[i, 'pdbSequence'] == 'nan':
            aligned_nm.at[i, 'mutationPositionOnPDB'] = 'nan'
            aligned_nm.at[i, 'domainStartonPDB'] = 'nan'
            aligned_nm.at[i, 'domainEndonPDB'] = 'nan'
            aligned_nm.at[i, 'pdb_alignStatus'] = 'nan'

    # Check if they the same column name before merging.
    aligned_m = aligned_m.astype(str)
    aligned_nm = aligned_nm.astype(str)

    # aligned_m.to_csv('/Users/fatmacankara/Desktop/new_benchmark/aligned_m.txt', sep='\t', index=False)
    # aligned_nm.to_csv('/Users/fatmacankara/Desktop/new_benchmark/aligned_nm.txt', sep='\t', index=False)

    frames = [aligned_m, aligned_nm]
    after_up_pdb_alignment = pd.concat(frames, sort=False)
    if len(after_up_pdb_alignment) == 0:
        after_up_pdb_alignment['pdb_alignStatus'] = ''
        after_up_pdb_alignment['mutationPositionOnPDB'] = ''
        after_up_pdb_alignment['domainStartonPDB'] = ''
        after_up_pdb_alignment['domainEndonPDB'] = ''

    after_up_pdb_alignment = after_up_pdb_alignment.sort_values(
        by=['uniprotID', 'wt', 'mut', 'pos', 'pdb_alignStatus', 'resolution', 'chain'],
        ascending=[True, True, True, True, True, True, True])

    after_up_pdb_alignment = after_up_pdb_alignment.drop_duplicates(['uniprotID', 'wt', 'mut', 'pos'], keep='first')

    after_up_pdb_alignment = after_up_pdb_alignment.astype('str')

    pdb_aligned = after_up_pdb_alignment[
        (after_up_pdb_alignment.pdbID != 'nan') & (after_up_pdb_alignment.mutationPositionOnPDB != 'nan')]
    yes_pdb_no_match = after_up_pdb_alignment[
        (after_up_pdb_alignment.pdbID != 'nan') & (after_up_pdb_alignment.mutationPositionOnPDB == 'nan')]
    no_pdb = no_pdb.copy()

    # len(aligned.drop_duplicates(['datapoint'])) + len(no_pdb.drop_duplicates(['datapoint'])) + len(yes_pdb_no_match.drop_duplicates(['datapoint'])) + len(not_match_in_uniprot.drop_duplicates(['datapoint'])) == len(data)

    print('PDB matching is completed...\n')
    print('SUMMARY')
    print('-------')
    print('%d data points that failed to match a UniProt Sequence are discarded.' % len(
        not_match_in_uniprot.drop_duplicates(['datapoint'])))
    print('Of the remaining %d:' % uniprot_matched_size)
    print('--%d of %d successfully aligned with PDB structures.' % (
        len(pdb_aligned.drop_duplicates(['datapoint'])), with_pdb_size))
    print('--%d of %d not found on the covered area by the structure.' % (
        len(yes_pdb_no_match.drop_duplicates(['datapoint'])), with_pdb_size))
    print('--PDB structures not found for %d datapoints.' % len(no_pdb.drop_duplicates(['datapoint'])))
    print('--%d will be searched in Swiss-Model database.\n' % (
            len(yes_pdb_no_match.drop_duplicates(['datapoint'])) + len(no_pdb.drop_duplicates(['datapoint']))))


    dfM = None
    dfNM = None
    aligned_nm = None
    aligned_m = None
    after_up_pdb_alignment = None

    print('Proceeding to  SwissModel search...')
    print('------------------------------------\n')

    # At this point we have 4 dataframes
    # 1. after_up_pdb_alignment --- This is after PDB sequence alignment. There may be mutations that wasnt found matching to after the alignment. Will be searched in other databases as well.
    # 1a. aligned --- we are done with this.
    # 1b. yes_pdb_no_match --- They have PDB structures but not matched, so will be searched in the other databases.
    # 2. not_match_in_uniprot --- This wont be aligned with anything because these proteins dont have a uniprot ID. Only basic info is present.
    # 3. no_pdb --- No PDB structures were found for them. Will be searched in other databases.

    """
    Step 8
    Neutralize data points that are to be searched in Swiss-Model
    # One point is that yes_pdb_no_match's annotations are the adjusted according to the PDBs they are matched before.
    # They need to be converted to their old original UniProt annotation positions.
    """
    yes_pdb_no_match.drop(['disulfide', 'intMet',
                           'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                           'nucleotideBinding', 'lipidation', 'site', 'transmembrane', 'crosslink',
                           'mutagenesis', 'strand', 'helix', 'turn', 'metalBinding', 'repeat',
                           'caBinding', 'topologicalDomain', 'bindingSite', 'region',
                           'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil',
                           'peptide', 'transitPeptide', 'glycosylation', 'propeptide', 'disulfideBinary',
                           'intMetBinary', 'intramembraneBinary',
                           'naturalVariantBinary', 'dnaBindingBinary', 'activeSiteBinary',
                           'nucleotideBindingBinary', 'lipidationBinary', 'siteBinary',
                           'transmembraneBinary', 'crosslinkBinary', 'mutagenesisBinary',
                           'strandBinary', 'helixBinary', 'turnBinary', 'metalBindingBinary',
                           'repeatBinary', 'topologicalDomainBinary', 'caBindingBinary',
                           'bindingSiteBinary', 'regionBinary', 'signalPeptideBinary',
                           'modifiedResidueBinary', 'zincFingerBinary', 'motifBinary',
                           'coiledCoilBinary', 'peptideBinary', 'transitPeptideBinary',
                           'glycosylationBinary', 'propeptideBinary', 'pdbSequence', 'pdbInfo', 'pdbID',
                           'chain', 'resolution', 'pdb_alignStatus', 'mutationPositionOnPDB',
                           'domainStartonPDB', 'domainEndonPDB'], axis=1, inplace=True)

    to_swiss = pd.concat([yes_pdb_no_match.drop_duplicates(['datapoint']), no_pdb.drop_duplicates(['datapoint'])])
    no_pdb = None
    to_swiss.reset_index(inplace=True)
    to_swiss.drop(['index'], axis=1, inplace=True)
    to_swiss = to_swiss.astype('str')
    to_swiss = to_swiss.replace({'NaN': 'nan'})
    # to_swiss.to_csv('/Users/fatmacankara/Desktop/new_benchmark/to_swiss.txt', sep='\t', index=False)

    # Modellerin özeti olan dosyayı oluştur.
    if len(to_swiss) != 0:
        print('Generating SwissModel file...\n')

        swiss_model = pd.read_csv(swiss_model_path, sep='\t',
                                  dtype=str, header=None, skiprows=1,
                                  names=['input', 'UniProtKB_ac', 'iso_id', 'uniprot_seq_length', 'uniprot_seq_md5',
                                         'coordinate_id', 'provider', 'from', 'to', 'template', 'qmean', 'qmean_norm',
                                         'seqid', 'url']
                                  )
        for i in swiss_model.index:
            k = ((swiss_model.at[i, 'input']).split(','))
            for j in k:
                if ('uniprot_ac' in str(j)):
                    swiss_model.at[i, 'UniProtKB_ac'] = j.split(':')[1].strip('\"')
                if 'iso_id' in str(j):
                    swiss_model.at[i, 'iso_id'] = j.split(':')[1].strip('\"')
                if 'uniprot_seq_length' in str(j):
                    swiss_model.at[i, 'uniprot_seq_length'] = j.split(':')[1].strip('\"')
                if 'uniprot_seq_md5' in str(j):
                    swiss_model.at[i, 'uniprot_seq_md5'] = j.split(':')[1].strip('\"')
                if 'coordinate_id' in str(j):
                    swiss_model.at[i, 'coordinate_id'] = j.split(':')[1].strip('\"')
                if 'provider' in str(j) and 'http' not in str(j):
                    swiss_model.at[i, 'provider'] = j.split(':')[1].strip('\"')
                if 'from' in str(j) and 'http' not in str(j):
                    swiss_model.at[i, 'from'] = j.split(':')[1].strip('\"')
                if 'to' in str(j) and 'http' not in str(j):
                    swiss_model.at[i, 'to'] = j.split(':')[1].strip('\"')
                if 'template' in str(j) and 'http' not in str(j):
                    swiss_model.at[i, 'template'] = j.split(':')[1].strip('\"')
                if 'qmean' in str(j) and 'norm' not in str(j):
                    swiss_model.at[i, 'qmean'] = j.split(':')[1].strip('\"')
                if 'qmean_norm' in str(j):
                    swiss_model.at[i, 'qmean_norm'] = j.split(':')[1].strip('\"')
                if 'seqid' in str(j):
                    swiss_model.at[i, 'seqid'] = j.split(':')[1].strip('\"')
                if 'url' in str(j):
                    swiss_model.at[i, 'url'] = ''.join(j.split(':')[1:])
    else:
        swiss_model = pd.DataFrame(
            columns=['input', 'UniProtKB_ac', 'iso_id', 'uniprot_seq_length', 'uniprot_seq_md5', 'coordinate_id',
                     'provider', 'from', 'to', 'template', 'qmean', 'qmean_norm', 'seqid', 'url', 'whichIsoform'])

    swiss_model = swiss_model.astype('str')
    try:
        swiss_model.iso_id = swiss_model.iso_id.astype('str')
    except:
        AttributeError
        swiss_model['iso_id'] = 'nan'
    swiss_model = swiss_model[swiss_model.UniProtKB_ac != 'nan']
    for ind in swiss_model.index:
        swiss_model.at[ind, 'UniProtKB_ac'] = swiss_model.at[ind, 'UniProtKB_ac'].split('-')[0]
        if swiss_model.at[ind, 'iso_id'] != 'nan':
            swiss_model.at[ind, 'whichIsoform'] = swiss_model.at[ind, 'iso_id'].split('-')[1]
        else:
            swiss_model.at[ind, 'whichIsoform'] = 'nan'
    swiss_model.drop(['input'], axis=1, inplace=True)
    swiss_model = swiss_model[swiss_model.provider == 'SWISSMODEL']
    print('Index File Processed...\n')
    # Get relevant columns
    swiss_model = swiss_model[['UniProtKB_ac', 'from', 'to', 'template', 'qmean_norm', 'seqid', 'url', 'whichIsoform']]
    # Sort models on qmean score and identity. Some proteins have more than one models, we will pick one.
    swiss_model = swiss_model.sort_values(by=['UniProtKB_ac', 'qmean_norm', 'seqid'], ascending=False)
    swiss_model.reset_index(inplace=True)
    swiss_model.drop(['index'], axis=1, inplace=True)

    # Get protein IDs for which there exist models.
    swiss_model_ids = set(swiss_model.UniProtKB_ac.to_list())
    to_swiss = to_swiss.astype(str)
    no_swiss_models = pd.DataFrame()
    for i in to_swiss.index:
        if to_swiss.at[i, 'uniprotID'] not in swiss_model_ids:
            k = pd.Series(to_swiss.iloc[i])
            no_swiss_models = no_swiss_models.append(k, ignore_index=True)

    no_swiss_models = no_swiss_models.astype(str)
    if len(no_swiss_models) == 0:
        no_swiss_models = pd.DataFrame(columns=to_swiss.columns)
    else:
        no_swiss_models = no_swiss_models[to_swiss.columns]
        no_swiss_models.reset_index(inplace=True)
        no_swiss_models.drop('index', axis=1, inplace=True)

    with_swiss_models = pd.concat([to_swiss, no_swiss_models]).drop_duplicates(['datapoint'], keep=False)
    with_swiss_models = with_swiss_models[to_swiss.columns]

    # Add model info.

    with_swiss_models = with_swiss_models.astype(str)
    swiss_model = swiss_model.astype(str)

    swiss_models_with_data = pd.merge(with_swiss_models, swiss_model, left_on=['uniprotID', 'whichIsoform'],
                                      right_on=['UniProtKB_ac', 'whichIsoform'],
                                      how='left')
    swiss_models_with_data = swiss_models_with_data.astype(str)
    swiss_models_with_data = swiss_models_with_data.sort_values(by=['uniprotID', 'wt', 'mut', 'pos', 'qmean_norm'],
                                                                ascending=False)
    swiss_models_with_data = swiss_models_with_data.drop_duplicates()
    swiss_models_with_data = swiss_models_with_data.drop(['UniProtKB_ac', 'seqid'], axis=1)
    swiss_models_with_data.pos = swiss_models_with_data.pos.astype('int')
    swiss_models_with_data = swiss_models_with_data.astype(str)

    # Get the ones in the list but without model url and add to the list to go to modbase.
    url_nan = swiss_models_with_data[swiss_models_with_data.url == 'nan']

    # Add this nan's to no_model. These will be searched in MODBASE because here they dont have urls.
    url_nan = url_nan.drop(['from', 'qmean_norm', 'template', 'to', 'url'], axis=1)

    no_swiss_models_2 = pd.concat([no_swiss_models, url_nan])
    swiss_models_with_data = swiss_models_with_data[swiss_models_with_data.url != 'nan']
    for i in swiss_models_with_data.index:
        try:
            swiss_models_with_data.at[i, 'chain'] = swiss_models_with_data.at[i, 'template'].split('.')[2]
            swiss_models_with_data.at[i, 'template'] = swiss_models_with_data.at[i, 'template'].split('.')[0]
        except:
            IndexError
    if len(swiss_models_with_data) == 0:
        swiss_models_with_data['chain'] = ''
        swiss_models_with_data['template'] = ''

    swiss_models_with_data.qmean_norm = swiss_models_with_data.qmean_norm.astype('str')
    swiss_models_with_data.chain = swiss_models_with_data.chain.astype('str')
    swiss_models_with_data['qmean_norm'] = swiss_models_with_data.qmean_norm.apply(lambda x: round(float(x), 2))
    swiss_models_with_data = swiss_models_with_data.astype(str)

    # swiss_models_with_data: These data points will be aligned with their corresponding model sequences.
    # Add sequences  # Burda sayı düşüyor.

    no_swiss_models_2.reset_index(inplace=True)
    no_swiss_models_2.drop('index', axis=1, inplace=True)

    swiss_models_with_data.reset_index(inplace=True)
    swiss_models_with_data.drop('index', axis=1, inplace=True)
    swiss_model_ids = None
    with_swiss_models = None
    swiss_model = None
    no_swiss_models = None
    url_nan = None
    # len(to_swiss.drop_duplicates(['datapoint'])) == len(no_swiss_models_2.drop_duplicates(['datapoint'])) + len(swiss_models_with_data.drop_duplicates(['datapoint']))

    # At this point we have:
    # pdb_aligned --- Align in the PDB phase
    # not_match_in_uniprot --- This wont be aligned with anything because these proteins dont have a uniprot ID. Only basic info is present.
    # to_swiss (no_pdb + yes_pdb_no_match) --- to be searched in SwissModel database
    # to_swiss (with_swiss_models & no_swiss_models)
    # swiss_models_with_data --- We found swiss models for them.
    # no_swiss_models_2 (no_swiss_models + url_nan)--- to be searched in modbase (the ones having swissmodels but not matching with the boundaries  & broken_swiss will be added here)

    """
    STEP 9
    Associated model IDs are added. 
    Download model files.
    """
    print('Beginning SwissModel files download...')
    existing_swiss = glob.glob(path_to_output_files + 'swissmodel_structures/*')
    existing_swiss = ['.'.join(i.split('/')[-1].split('.')[:-1]) for i in existing_swiss]
    swissmodels_fasta = pd.DataFrame()
    for i in swiss_models_with_data.index:
        protein = swiss_models_with_data.at[i, 'uniprotID']
        template = swiss_models_with_data.at[i, 'template'].split('.')[0]
        qmean_norm = str(round(float(swiss_models_with_data.at[i, 'qmean_norm']), 2))
        if protein + '_' + template + '_' + qmean_norm not in existing_swiss:
            url = swiss_models_with_data.at[i, 'url'].strip('\"').strip('}').replace('\\', '').strip('\"').replace(
                'https',
                'https:')
            req = r.get(url)
            name = path_to_output_files + 'swissmodel_structures/' + protein + '_' + template + '_' + qmean_norm + '.txt'
            print('Downloading for Protein:', protein + ' Model: ' + template)
            with open(name, 'wb') as f:
                f.write(req.content)
        else:
            print('Model exists.')
            name = glob.glob(
                path_to_output_files + 'swissmodel_structures/' + protein + '_' + template + '_' + qmean_norm + '.txt')[
                0]
        with open(name, encoding="utf8") as f:
            fasta = ''
            lines = f.readlines()
            chain = ''
            for row in lines:
                if row[0:4] == 'ATOM' and row[13:15] == 'CA':
                    chain = row[20:22].strip()
                    fasta += threeToOne(row[17:20])
                if row[0:3] == 'TER':
                    k = pd.Series([protein, template, qmean_norm, chain.upper(), fasta])
                    swissmodels_fasta = swissmodels_fasta.append(k, ignore_index=True)
                    fasta = ''

    if len(swissmodels_fasta) == 0:
        swissmodels_fasta = pd.DataFrame(columns=['uniprotID', 'template', 'qmean_norm', 'chain', 'fasta'])
    else:
        # swissmodels_fasta = no_swiss_models[to_swiss.columns]  bu satır yanlış olabilir mi? bi anlamsız geldii. belki yer değiştirirsekn hata yaptım.25 mart.
        swissmodels_fasta.columns = ['uniprotID', 'template', 'qmean_norm', 'chain', 'fasta']
    swissmodels_fasta = swissmodels_fasta.astype(str)
    # swissmodels_fasta.to_csv('/Users/fatmacankara/Desktop/new_benchmark/swiss_models_fasta.txt', sep='\t', index=False)

    swiss_models_with_data.qmean_norm = swiss_models_with_data.qmean_norm.astype(float)
    swissmodels_fasta.qmean_norm = swissmodels_fasta.qmean_norm.astype(float)

    # Bu dosya üerince şunu yap; bazı aynı modellerin sekansı yok.
    swissmodels_fasta = swissmodels_fasta.sort_values(['uniprotID', 'template', 'qmean_norm', 'chain'],
                                                      axis=0)  # example = 3gdh
    swissmodels_fasta.reset_index(inplace=True)
    swissmodels_fasta.drop(['index'], axis=1, inplace=True)
    swissmodels_fasta = swissmodels_fasta.drop_duplicates(['uniprotID', 'template', 'qmean_norm', 'chain'])
    swissmodels_fasta = swissmodels_fasta.drop_duplicates(['uniprotID', 'template', 'chain', 'fasta'])
    swissmodels_fasta = swissmodels_fasta.drop_duplicates(['uniprotID', 'template', 'fasta'])
    # Some files were broken, thus their PDBs couldnt be recorded.
    swissmodels_fasta = swissmodels_fasta.drop_duplicates()
    swissmodels_fasta = swissmodels_fasta.astype(str)

    swiss_models_with_data = swiss_models_with_data.astype(str)
    swissmodels_fasta = swissmodels_fasta.astype(str)
    swiss_models_with_data1 = swiss_models_with_data.merge(swissmodels_fasta,
                                                           on=['uniprotID', 'template', 'qmean_norm', 'chain'])
    swiss_models_with_data1 = swiss_models_with_data1.sort_values(['datapoint', 'fasta'], axis=0,
                                                                  ascending=[True, False])
    swiss_models_with_data1 = swiss_models_with_data1.drop_duplicates(['datapoint', 'template'])

    # Bazı proteinlere model gelmedi, bunları da modbasee gideceklere ekle. Bunu oluşturmak çok uzun sürdü. O yüzden yazdırıyorum.
    # swiss_models_with_data1'de sadece tam model bilgisi olanlar var.
    swiss_models_with_data1_dp = list(set(swiss_models_with_data1.datapoint.to_list()))
    swiss_models_with_data.reset_index(inplace=True)
    swiss_models_with_data.drop(['index'], axis=1, inplace=True)

    broken_swiss = pd.DataFrame()
    c = 0
    for i in swiss_models_with_data.index:  # en baştaki dfde var ama model gelende yok.
        if swiss_models_with_data.at[i, 'datapoint'] not in swiss_models_with_data1_dp:
            k = pd.Series(swiss_models_with_data.iloc[i])
            broken_swiss = broken_swiss.append(k, ignore_index=True)
            c += 1
    if len(broken_swiss) == 0:
        broken_swiss = pd.DataFrame(columns=swiss_models_with_data.columns.to_list())

    swiss_models_with_data = swiss_models_with_data1.copy()
    # Bundan sonra bile bazıları hem sekansı olan hem de sekansı olmayan IDlerle eşlenmiş. Bu durumda sekansın olmadığı durumları yok etmek lazım.
    # Çünkü eşlenmezse zaten sonraki adıma alıcaz.
    swiss_models_with_data.qmean_norm = swiss_models_with_data.qmean_norm.astype('float')
    swiss_models_with_data = swiss_models_with_data.sort_values(['uniprotID', 'wt', 'mut', 'qmean_norm'],
                                                                axis=0, ascending=[True, True, True, False])

    # Delete the same model sequence with lower quality
    swiss_models_with_data = swiss_models_with_data.drop_duplicates(['uniprotID', 'wt', 'mut', 'pos', 'fasta'],
                                                                    keep='first')
    swiss_models_with_data.uniprotSequence = swiss_models_with_data.uniprotSequence.astype('str')
    swiss_models_with_data.pos = swiss_models_with_data.pos.astype('int')
    len(swiss_models_with_data.drop_duplicates(['datapoint'])) + len(broken_swiss.drop_duplicates(['datapoint'])) + len(
        no_swiss_models_2.drop_duplicates(['datapoint'])) == len(to_swiss.drop_duplicates(['datapoint']))
    # This printed data here includes all possible models with different qualities,
    # because we may get a hit in either of them.
    swiss_models_with_data.rename({'fasta': 'pdbSequence'}, axis=1, inplace=True)  # for convenience.

    # NOW DO ALIGNMENT HERE

    swiss_models_with_data = swiss_models_with_data.replace({'[\'?\']': 'nan'})
    swiss_models_with_data = swiss_models_with_data.replace({'[]': 'nan'})
    swiss_models_with_data.rename({'template': 'pdbID'}, axis=1,
                                  inplace=True)  # Only to be able use the alignment code above.
    swiss_models_with_data = swiss_models_with_data.astype(str)
    swiss_models_with_data.pdbSequence = swiss_models_with_data.pdbSequence.astype('str')
    swiss_models_with_data = add_annotations(swiss_models_with_data)
    swiss_models_with_data = swiss_models_with_data.astype(str)
    swiss_models_with_data.replace({'NaN': 'nan'}, inplace=True)
    # swiss_models_with_data.to_csv('/Users/fatmacankara/Desktop/new_benchmark/right_before_alignment.txt', sep='\t', index=False)
    swiss_models_with_data_copy = swiss_models_with_data.copy()
    #  burda alignment zaten ana dataframei de döndürüyo o yüzden kopya al bunlan önce ve onunla karşılaştır değişimi.
    swiss_models_with_data1_dp = None
    swiss_models_with_data1 = None
    existing_swiss = None
    swissmodels_fasta = None

    print('Aligning sequences...\n')
    swiss_models_with_data['uniprotSequence'] = swiss_models_with_data['uniprotSequence'].str.replace('U', 'C')
    swiss_models_with_data['pdbSequence'] = swiss_models_with_data['pdbSequence'].str.replace('U', 'C')
    swiss_model_aligned = alignment(swiss_models_with_data)
    swiss_models_with_data = None

    # swiss_model_aligned.to_csv('/Users/fatmacankara/Desktop/new_benchmark/swiss_alignment.txt', sep='\t', index=False)

    if len(swiss_model_aligned) == 0:
        swiss_model_aligned = pd.DataFrame(columns=pdb_aligned.columns)
        swiss_model_aligned['qmean_norm'] = 'nan'
    else:
        swiss_model_aligned = swiss_model_aligned.astype(str)
        swiss_model_aligned.replace({'NaN': 'nan'}, inplace=True)


    # Some datapoints appear in both nan and not_nan. If not_nan we take it only once.
    nan = swiss_model_aligned[swiss_model_aligned.mutationPositionOnPDB == 'nan']
    not_nan = swiss_model_aligned[swiss_model_aligned.mutationPositionOnPDB != 'nan']
    not_nan.qmean_norm = not_nan.qmean_norm.astype('float')
    not_nan.sort_values(['datapoint', 'pdb_alignStatus', 'qmean_norm'], ascending=[True, True, False], inplace=True)

    which_ones_are_match = pd.concat([not_nan, nan]).drop_duplicates(['datapoint'], keep='first')
    swiss_match = which_ones_are_match[which_ones_are_match.mutationPositionOnPDB != 'nan']
    swiss_not_match = which_ones_are_match[which_ones_are_match.mutationPositionOnPDB == 'nan']

    swiss_match.qmean_norm = swiss_match.qmean_norm.astype('float')
    swiss_match.sort_values(['uniprotID', 'wt', 'pos', 'mut', 'pdb_alignStatus', 'qmean_norm'],
                            ascending=[True, True, True, True, True, False], inplace=True)
    swiss_match.drop_duplicates(['uniprotID', 'wt', 'pos', 'mut'], keep='first', inplace=True)

    # swiss_match.to_csv('/Users/fatmacankara/Desktop/benchmark_dec2020/swiss_aligned_match.txt', sep=' ', index=False)

    # swiss_not_match.to_csv('/Users/fatmacankara/Desktop/benchmark_dec2020/swiss_aligned_not_match_to_modbase_ekle.txt', sep=' ', index=False)
    swiss_not_match = swiss_not_match[no_swiss_models_2.columns]
    broken_swiss = broken_swiss[no_swiss_models_2.columns]
    swiss_not_match = swiss_not_match.drop_duplicates(['datapoint'])
    broken_swiss = broken_swiss.drop_duplicates(['datapoint'])

    to_modbase = pd.concat([no_swiss_models_2, broken_swiss]).drop_duplicates()
    to_modbase = pd.concat([to_modbase, swiss_not_match]).drop_duplicates()
    to_modbase = to_modbase.astype(str)
    to_swiss_columns = to_swiss.columns
    to_swiss_size = len(to_swiss.drop_duplicates(['datapoint']))
    to_swiss = None
    # to_modbase.to_csv('/Users/fatmacankara/Desktop/new_benchmark/to_modbase.txt', sep = '\t', index=False)

    # CONTROL

    """
    # This should be the whole data.
    len(swiss_match.drop_duplicates(['datapoint'])) + len(aligned.drop_duplicates(['datapoint'])) + len(to_modbase.drop_duplicates(['datapoint'])) + len(not_match_in_uniprot.drop_duplicates(['datapoint'])) ,len(data)
    len(aligned.drop_duplicates(['datapoint'])) + len(not_match_in_uniprot.drop_duplicates(['datapoint'])) +len(to_swiss.drop_duplicates(['datapoint']))== len(data)
    """
    print('SwissModel matching is completed...\n')
    print('SUMMARY')
    print('-------')
    print('%d data points that failed to match a UniProt Sequence are discarded.' % len(
        not_match_in_uniprot.drop_duplicates(['datapoint'])))
    print('Of the remaining %d:' % uniprot_matched_size)
    print('--%d of %d successfully aligned with PDB structures.' % (
        len(pdb_aligned.drop_duplicates(['datapoint'])), with_pdb_size))
    print('--%d of %d successfully aligned with SwissModels structures.' % (
        len(swiss_match.drop_duplicates(['datapoint'])), to_swiss_size))
    print('--%d will be searched in ModBase database.\n' % len(to_modbase.drop_duplicates(['datapoint'])))

    print('Proceeding to ModBase search...')
    print('------------------------------------\n')
    no_swiss_models_2 = None
    broken_swiss = None
    swiss_model_aligned = None
    nan = None
    not_nan = None
    which_ones_are_match = None
    swiss_not_match = None

    # STEP :  GO TO MODBASE
    # Should not include anything related to prev models.
    if len(to_modbase) != 0:
        to_modbase = to_modbase.astype(str)

        # GET MODBASE MODELS

        # Get IDs from data to retrieve only their models from MODBASE
        to_modbase.reset_index(inplace=True)
        to_modbase.drop(['index'], axis=1, inplace=True)

        existing_modbase_models = glob.glob(path_to_output_files + 'modbase_structures/*')
        existing_modbase_models = [i.split('/')[-1].split('.')[0] for i in existing_modbase_models]
        existing_modbase_models_ind = glob.glob(path_to_output_files + 'modbase_structures_individual/*')
        existing_modbase_models_ind = [i.split('/')[-1].split('.')[0] for i in existing_modbase_models_ind]
        modbase_reduced = pd.DataFrame()
        modbase_fasta = pd.DataFrame()

        print('Retrieving ModBase models...\n')
        # Get model files associated with each UniProtID
        for protein in list(set(to_modbase.uniprotID.to_list())):
            if protein not in existing_modbase_models:
                print('Downloading Modbase models for ', protein)
                url = 'https://salilab.org/modbase/retrieve/modbase/?databaseID=' + protein
                req = r.get(url)
                name = path_to_output_files + 'modbase_structures/' + protein + '.txt'
                with open(name, 'wb') as f:
                    f.write(req.content)
            else:
                print('Model exists for', protein)
                name = glob.glob(path_to_output_files + 'modbase_structures/' + protein + '.txt')[0]
            with open(name, encoding="utf8") as f:
                a = open(name, 'r').read()
                soup = BeautifulSoup(a, 'lxml')
                for pdb in soup.findAll('pdbfile'):
                    model_id = str(pdb.contents[1])[10:-11]
                    if model_id not in existing_modbase_models_ind:
                        with open(path_to_output_files + 'modbase_structures_individual/' + model_id + '.txt', 'w',
                                  encoding="utf8") as individual:
                            individual.write(str('UniProt ID: ' + protein))
                            individual.write('\n')
                            individual.write(str(pdb.contents[3])[10:-11].strip())
                    with open(path_to_output_files + 'modbase_structures_individual/' + model_id + '.txt',
                              encoding="utf8") as f:
                        fasta = ''
                        chain = ''
                        template_chain = ''
                        score = -999
                        for ind_line in f.readlines():
                            if ind_line[0:10] == 'UniProt ID':
                                uniprot_id = ind_line.split(':')[1].strip()
                            if ind_line[0:23] == 'REMARK 220 TARGET BEGIN':
                                target_begin = ind_line[40:43].strip()
                            if ind_line[0:21] == 'REMARK 220 TARGET END':
                                target_end = ind_line[40:43].strip()
                            if ind_line[0:25] == 'REMARK 220 TEMPLATE BEGIN':
                                pdb_begin = ind_line[40:43].strip()
                            if ind_line[0:23] == 'REMARK 220 TEMPLATE END':
                                pdb_end = ind_line[40:43].strip()
                            if ind_line[0:23] == 'REMARK 220 TEMPLATE PDB':
                                pdb_code = ind_line[40:43].strip()
                            if ind_line[0:25] == 'REMARK 220 TEMPLATE CHAIN':
                                pdb_chain = ind_line[40:43].strip()
                            if ind_line[0:32] == 'REMARK 220 ModPipe Quality Score':
                                quality_score = ind_line[40:].strip()
                            if ind_line[0:27] == 'REMARK 220 MODPIPE MODEL ID':
                                model_id = ind_line[40:].strip()
                            if ind_line[0:25] == 'REMARK 220 TEMPLATE CHAIN':
                                template_chain = ind_line[40:42].strip()
                            if ind_line[0:4] == 'ATOM' and ind_line[13:15] == 'CA':
                                fasta += threeToOne(ind_line[17:20])
                            if ind_line[0:32] == 'REMARK 220 ModPipe Quality Score':
                                try:
                                    score = ind_line[40:].strip()
                                except (ValueError):
                                    score = -999
                            if ind_line[0:3] == 'TER' or ind_line[0:3] == 'END':
                                k = pd.Series([uniprot_id, model_id, str(score), template_chain, fasta])
                                modbase_fasta = modbase_fasta.append(k, ignore_index=True)
                                fasta = ''
                        try:
                            k = pd.Series(
                                [uniprot_id, target_begin, target_end, pdb_code, pdb_chain, pdb_begin, pdb_end,
                                 quality_score,
                                 model_id])
                            modbase_reduced = modbase_reduced.append(k, ignore_index=True)
                        except:
                            NameError
                            print('This file doesnt have Quality Score. Replacer: -999', model_id)
                            quality_score = -999

        print()
        if len(modbase_fasta) != 0:
            modbase_fasta.columns = ['uniprotID', 'template', 'score', 'chain', 'fasta']
        else:
            modbase_fasta = pd.DataFrame(columns=['uniprotID', 'template', 'score', 'chain', 'fasta'])
        modbase_fasta = modbase_fasta.astype(str)
        modbase_fasta = modbase_fasta.replace({'': 'nan'})
        modbase_fasta = modbase_fasta.replace({'NaN': 'nan'})
        modbase_fasta = modbase_fasta[modbase_fasta.fasta != 'nan']

        print('Modbase model frame constructed.\n')
        if len(modbase_reduced) != 0:
            modbase_reduced.columns = ['UniprotID', 'TargetBeg', 'TargetEnd', 'PDBCode', 'PDBChain', 'PDBBegin',
                                       'PDBEnd',
                                       'ModPipeQualityScore', 'ModelID']
        else:
            modbase_reduced = pd.DataFrame(
                columns=['UniprotID', 'TargetBeg', 'TargetEnd', 'PDBCode', 'PDBChain', 'PDBBegin', 'PDBEnd',
                         'ModPipeQualityScore', 'ModelID'])

        # modbase_reduced.to_csv('/Users/fatmacankara/Desktop/new_benchmark/modbase_reduced.txt', index=False, sep=' ')
        to_modbase = add_annotations(to_modbase)
        to_modbase = to_modbase.astype(str)
        to_modbase.fillna('nan', inplace=True)
        to_modbase = to_modbase.replace({'NaN': 'nan'})
        to_modbase.replace({'[]': 'nan'}, inplace=True)
        to_modbase.replace({'nan-nan': 'nan'}, inplace=True)
        to_modbase.replace({'': 'nan'}, inplace=True)

        # to_modbase.to_csv(
        #    '/Users/fatmacankara/Desktop/varibench_intermediate_files/to_modbase_varibench_annotations_added.txt', sep='\t',
        #    index=False)
        model_info_added = to_modbase.merge(modbase_reduced, right_on='UniprotID', left_on='uniprotID',
                                            how='left')
        modbase_reduced = None
        existing_modbase_models = None
        existing_modbase_models_ind = None


        model_info_added = model_info_added.drop(['UniprotID'], axis=1)
        model_info_added = model_info_added.rename(columns={'TargetBeg': 'from', 'TargetEnd': 'to',
                                                            'PDBCode': 'template', 'PDBChain': 'chain',
                                                            'ModPipeQualityScore': 'score',
                                                            'ModelID': 'pdbID'})
        model_info_added.drop(['PDBEnd', 'PDBBegin'], axis=1, inplace=True)
        model_info_added.score = model_info_added.score.astype(float)
        model_info_added = model_info_added.sort_values(by=['datapoint', 'score'],
                                                        ascending=False)
        model_info_added.reset_index(inplace=True)
        model_info_added.drop(['index'], axis=1, inplace=True)
        model_info_added = model_info_added.drop_duplicates()

        model_info_added = model_info_added.astype(str)
        # model_info_added.to_csv('/Users/fatmacankara/Desktop/new_benchmark/modbaseInfoAdded.txt', sep='\t', index=False)
        model_info_added = model_info_added.replace({'NaN': 'nan'})
        no_info = model_info_added[model_info_added.pdbID == 'nan']  # Model gelmedi, inmedi.
        with_modbase_info = model_info_added[model_info_added.pdbID != 'nan']  # Model gelenler.
        model_info_added = None

        len(no_info.drop_duplicates(['datapoint'])), len(with_modbase_info.drop_duplicates(['datapoint']))
        len(no_info.drop_duplicates(['datapoint'])) + len(with_modbase_info.drop_duplicates(['datapoint'])) == len(
            to_modbase.drop_duplicates(['datapoint']))

        # Add no_info to the rest down below!
        no_info = no_info[to_swiss_columns]

        with_modbase_info.score = with_modbase_info.score.astype(float)
        modbase_fasta.score = modbase_fasta.score.astype(float)

        # Bu dosya üerince şunu yap; bazı aynı modellerin sekansı yok.
        modbase_fasta = modbase_fasta.sort_values(['uniprotID', 'score', 'template', 'chain'],
                                                  ascending=[True, False, True, True], axis=0)  # example = 3gdh
        # modbase_fasta.to_csv('/Users/fatmacankara/Desktop/new_benchmark/modbase_fasta.txt', sep='\t', index=False)

        # I added this newly downloaded ones to the main model file.

        # Gereksiz bir ekstra column ama olmazsa yanlış okıuyo, yani sil.
        modbase_fasta = modbase_fasta.rename(columns={'template': 'pdbID'})
        with_modbase_info.pos = with_modbase_info.pos.astype('int')
        with_modbase_info.score = with_modbase_info.score.astype(float)
        with_modbase_info.score = with_modbase_info.score.apply(lambda x: round(x, 2))
        modbase_fasta.score = modbase_fasta.score.astype(float)
        modbase_fasta.score = modbase_fasta.score.apply(lambda x: round(x, 2))

        with_modbase_info = with_modbase_info.merge(modbase_fasta, on='pdbID', how='left')

        with_modbase_info.drop(['score_y'], axis=1, inplace=True)
        with_modbase_info.rename(columns={'score_x': 'score'}, inplace=True)
        with_modbase_info.drop(['uniprotID_y', 'chain_y'], axis=1, inplace=True)
        with_modbase_info.rename(columns={'uniprotID_x': 'uniprotID', 'chain_x': 'chain'}, inplace=True)

        with_modbase_info.score = with_modbase_info.score.astype('float')
        with_modbase_info = with_modbase_info.sort_values(['uniprotID', 'wt', 'mut', 'pos', 'score', 'from', 'to'],
                                                          axis=0,
                                                          ascending=[True, True, True, True, False, True, False])
        with_modbase_info = with_modbase_info.drop_duplicates(['uniprotID', 'wt', 'mut', 'pos', 'fasta'], keep='first')

        with_modbase_info = with_modbase_info.replace({'[\'?\']': 'nan'})
        with_modbase_info = with_modbase_info.replace({'[]': 'nan'})
        with_modbase_info = with_modbase_info.replace({'\'?\', ': ''})
        with_modbase_info = with_modbase_info.replace({', \'?\'': ''})
        with_modbase_info = with_modbase_info.replace({'(': ''})
        with_modbase_info = with_modbase_info.replace(
            {')': ''})  # Bu genelde olmuyo neden bilmiyorum. Manual olarak düzelttim.
        # Ve düzeltilmiş dosyayı asıl dosyaymış gibi kaydettim. Bu muhtemelen anotasyonları eklediğimiz yerle alakalı.
        with_modbase_info = with_modbase_info.astype(str)
        with_modbase_info.fasta = with_modbase_info.fasta.astype('str')
        with_modbase_info.reset_index(inplace=True)
        with_modbase_info.drop('index', axis=1, inplace=True)


        align = with_modbase_info[
            with_modbase_info.fasta != 'nan']  # burda bir proteinden bir sürü var sıralanmış halde.
        yes_pdb_no_match = with_modbase_info[
            with_modbase_info.fasta == 'nan']  # modeli geldi ama gelen modelin sekansı yok. burda bir protein için x modelnin fastası olmayabilitr ama aynı protein diğer dataframede de görülebilir çünkü model y için fasta gelmiştir mesela.
        yes_pdb_no_match = yes_pdb_no_match[~yes_pdb_no_match.datapoint.isin(align.datapoint.to_list())]

        align.rename(columns={'fasta': 'pdbSequence'}, inplace=True)
        align['uniprotSequence'] = align['uniprotSequence'].str.replace('U', 'C')
        align['pdbSequence'] = align['pdbSequence'].str.replace('U', 'C')

        to_modbase_size = len(to_modbase.drop_duplicates(['datapoint']))
        modbase_fasta = None
        to_modbase = None

        print('Aligning sequences...\n')
        modbase_aligned = alignment(align)
        modbase_aligned = modbase_aligned.astype(str)
        modbase_aligned = modbase_aligned.replace({'NaN': 'nan'})
        # modbase_aligned.to_csv('/Users/fatmacankara/Desktop/new_benchmark/modbase_aligned.txt', sep='\t', index=False)

        # not_present_in_aligned = with_modbase_info[~with_modbase_info.datapoint.isin(modbase_aligned.datapoint.to_list())]

        # Get the ones whose models couldn't be found. Add to no_modbase (yani hiçbir şey de eşleşmemiş artık.)
        if len(with_modbase_info) != 0:
            not_in_aligned = pd.concat([modbase_aligned.drop_duplicates(['datapoint']),
                                        with_modbase_info.drop_duplicates(['datapoint'])]).drop_duplicates(
                ['datapoint'],
                keep=False)
        else:
            not_in_aligned = pd.DataFrame(columns=['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                                                   'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence',
                                                   'wt_sequence_match', 'whichIsoform', 'datapoint', 'disulfide',
                                                   'intMet',
                                                   'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                                                   'nucleotideBinding', 'lipidation', 'site', 'transmembrane',
                                                   'crosslink',
                                                   'mutagenesis', 'strand', 'helix', 'turn', 'metalBinding', 'repeat',
                                                   'topologicalDomain', 'caBinding', 'bindingSite', 'region',
                                                   'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif',
                                                   'coiledCoil',
                                                   'peptide', 'transitPeptide', 'glycosylation', 'propeptide',
                                                   'disulfide',
                                                   'intMet', 'intramembrane', 'naturalVariant', 'dnaBinding',
                                                   'activeSite',
                                                   'nucleotideBinding', 'lipidation', 'site', 'transmembrane',
                                                   'crosslink',
                                                   'mutagenesis', 'strand', 'helix', 'turn', 'metalBinding', 'repeat',
                                                   'topologicalDomain', 'caBinding', 'bindingSite', 'region',
                                                   'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif',
                                                   'coiledCoil',
                                                   'peptide', 'transitPeptide', 'glycosylation', 'propeptide', 'from',
                                                   'to', 'template', 'chain', 'score', 'pdbID', 'pdbSequence', 'fasta'])
        with_modbase_info = None
        if len(not_in_aligned) != 0:
            not_models = pd.concat([yes_pdb_no_match.drop_duplicates(['datapoint']),
                                    not_in_aligned.drop_duplicates(['datapoint'])]).drop_duplicates(['datapoint'],
                                                                                                    keep='first')
        # Retain the best model among the aligned ones. # Burayı swissteki gibi düzelt. Her ikisinde de olan ID olabilir.
        else:
            not_models = pd.DataFrame(columns=not_in_aligned.columns)

        yes_pdb_no_match = None
        # # Some datapoints appear in both nan and not_nan. If not_nan we take it only once.
        modbase_aligned = modbase_aligned.astype(str)
        if len(modbase_aligned) != 0:
            nan = modbase_aligned[modbase_aligned.mutationPositionOnPDB == 'nan']
            not_nan = modbase_aligned[modbase_aligned.mutationPositionOnPDB != 'nan']
            not_nan.score = not_nan.score.astype(float)
            not_nan.sort_values(['datapoint', 'pdb_alignStatus', 'score'], ascending=[True, True, False], inplace=True)

            not_nan = not_nan.sort_values(['datapoint', 'mutationPositionOnPDB', 'score'],
                                          ascending=[True, True, False])
            not_nan = not_nan.drop_duplicates(['datapoint'], keep='first')
        else:
            nan = pd.DataFrame(columns=modbase_aligned.columns)
            not_nan = pd.DataFrame(columns=modbase_aligned.columns)
        modbase_aligned = None
        which_ones_are_match = pd.concat([not_nan, nan]).drop_duplicates(['datapoint'], keep='first')
        if len(which_ones_are_match) == 0:
            which_ones_are_match = pd.DataFrame(
                columns=['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                         'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence',
                         'wt_sequence_match', 'whichIsoform', 'datapoint', 'disulfide', 'intMet',
                         'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                         'nucleotideBinding', 'lipidation', 'site', 'transmembrane', 'crosslink',
                         'mutagenesis', 'strand', 'helix', 'turn', 'metalBinding', 'repeat',
                         'topologicalDomain', 'caBinding', 'bindingSite', 'region',
                         'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil',
                         'peptide', 'transitPeptide', 'glycosylation', 'propeptide',
                         'disulfideBinary', 'intMetBinary', 'intramembraneBinary',
                         'naturalVariantBinary', 'dnaBindingBinary', 'activeSiteBinary',
                         'nucleotideBindingBinary', 'lipidationBinary', 'siteBinary',
                         'transmembraneBinary', 'crosslinkBinary', 'mutagenesisBinary',
                         'strandBinary', 'helixBinary', 'turnBinary', 'metalBindingBinary',
                         'repeatBinary', 'topologicalDomainBinary', 'caBindingBinary',
                         'bindingSiteBinary', 'regionBinary', 'signalPeptideBinary',
                         'modifiedResidueBinary', 'zincFingerBinary', 'motifBinary',
                         'coiledCoilBinary', 'peptideBinary', 'transitPeptideBinary',
                         'glycosylationBinary', 'propeptideBinary', 'from', 'to', 'template',
                         'chain', 'score', 'pdbID', 'pdbSequence', 'pdb_alignStatus',
                         'mutationPositionOnPDB', 'domainStartonPDB', 'domainEndonPDB'])
            modbase_match = which_ones_are_match[which_ones_are_match.mutationPositionOnPDB != 'nan']
            modbase_not_match = which_ones_are_match[which_ones_are_match.mutationPositionOnPDB == 'nan']

        else:
            modbase_match = which_ones_are_match[which_ones_are_match.mutationPositionOnPDB != 'nan']
            modbase_not_match = which_ones_are_match[which_ones_are_match.mutationPositionOnPDB == 'nan']

        which_ones_are_match = None
        modbase_match.score = modbase_match.score.astype('float')
        modbase_match = modbase_match.sort_values(['datapoint', 'mutationPositionOnPDB', 'score'],
                                                  ascending=[True, True, False])
        modbase_match.drop_duplicates(['datapoint'], keep='first', inplace=True)
        not_nan = None
        nan = None

        # modbase_match.to_csv('/Users/fatmacankara/Desktop/new_benchmark/modbase_aligned_match.txt', sep='\t', index=False)

        # merge not_in_align and modbase_not_match as they were both excluded from modbase match.

        # EN baştan modbase modeli yok diye ayırdıklarımız
        no_info = no_info[to_swiss_columns]
        no_info = no_info.drop_duplicates()

        # Modbase'de model olup sekans olmayanlar:
        not_models = not_models[to_swiss_columns]
        not_models = not_models.drop_duplicates()

        # Modbase modeli ve sekansı olup, pdbde eşleşmeyenler
        modbase_not_match = modbase_not_match[to_swiss_columns]
        modbase_not_match = modbase_not_match.drop_duplicates()
        if len(not_in_aligned) != 0 and len(modbase_not_match) != 0 and len(no_info) != 0:
            rest = pd.concat([not_in_aligned, modbase_not_match, no_info])
        elif len(not_in_aligned) != 0 and len(modbase_not_match) != 0 and len(no_info) == 0:
            rest = pd.concat([not_in_aligned, modbase_not_match])
        elif len(not_in_aligned) == 0 and len(modbase_not_match) != 0 and len(no_info) != 0:
            rest = pd.concat([modbase_not_match, no_info])
        elif len(not_in_aligned) != 0 and len(modbase_not_match) == 0 and len(no_info) != 0:
            rest = pd.concat([not_in_aligned, no_info])
        elif len(not_in_aligned) != 0 and len(modbase_not_match) == 0 and len(no_info) == 0:
            rest = not_in_aligned
        elif len(not_in_aligned) == 0 and len(modbase_not_match) != 0 and len(no_info) == 0:
            rest = modbase_not_match
        elif len(not_in_aligned) == 0 and len(modbase_not_match) == 0 and len(no_info) != 0:
            rest = no_info
        else:
            rest = pd.DataFrame(columns=['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                                         'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence',
                                         'wt_sequence_match', 'whichIsoform', 'datapoint'])

        rest = rest[to_swiss_columns]
        rest = rest.drop_duplicates()

        rest.reset_index(inplace=True)
        rest.drop(['index'], axis=1, inplace=True)
        rest = rest.astype('str')
        # est = add_annotations(rest)
        # pdbyle alakalı cols yok, concat yapucaksan ekle nan olarak.
        # rest.to_csv('/Users/fatmacankara/Desktop/new_benchmark/rest.txt', sep='\t', index=False)

    else:

        modbase_match = pd.DataFrame(columns=['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                                              'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence',
                                              'wt_sequence_match', 'whichIsoform', 'datapoint', 'disulfide', 'intMet',
                                              'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                                              'nucleotideBinding', 'lipidation', 'site', 'transmembrane', 'crosslink',
                                              'mutagenesis', 'strand', 'helix', 'turn', 'metalBinding', 'repeat',
                                              'topologicalDomain', 'caBinding', 'bindingSite', 'region',
                                              'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil',
                                              'peptide', 'transitPeptide', 'glycosylation', 'propeptide',
                                              'disulfideBinary', 'intMetBinary', 'intramembraneBinary',
                                              'naturalVariantBinary', 'dnaBindingBinary', 'activeSiteBinary',
                                              'nucleotideBindingBinary', 'lipidationBinary', 'siteBinary',
                                              'transmembraneBinary', 'crosslinkBinary', 'mutagenesisBinary',
                                              'strandBinary', 'helixBinary', 'turnBinary', 'metalBindingBinary',
                                              'repeatBinary', 'topologicalDomainBinary', 'caBindingBinary',
                                              'bindingSiteBinary', 'regionBinary', 'signalPeptideBinary',
                                              'modifiedResidueBinary', 'zincFingerBinary', 'motifBinary',
                                              'coiledCoilBinary', 'peptideBinary', 'transitPeptideBinary',
                                              'glycosylationBinary', 'propeptideBinary', 'from', 'to', 'template',
                                              'chain', 'score', 'pdbID', 'pdbSequence', 'pdb_alignStatus',
                                              'mutationPositionOnPDB', 'domainStartonPDB', 'domainEndonPDB'])
        not_in_aligned = pd.DataFrame(columns=['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                                               'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence',
                                               'wt_sequence_match', 'whichIsoform', 'datapoint', 'disulfide', 'intMet',
                                               'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                                               'nucleotideBinding', 'lipidation', 'site', 'transmembrane', 'crosslink',
                                               'mutagenesis', 'strand', 'helix', 'turn', 'metalBinding', 'repeat',
                                               'topologicalDomain', 'caBinding', 'bindingSite', 'region',
                                               'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil',
                                               'peptide', 'transitPeptide', 'glycosylation', 'propeptide', 'disulfide',
                                               'intMet', 'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                                               'nucleotideBinding', 'lipidation', 'site', 'transmembrane', 'crosslink',
                                               'mutagenesis', 'strand', 'helix', 'turn', 'metalBinding', 'repeat',
                                               'topologicalDomain', 'caBinding', 'bindingSite', 'region',
                                               'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil',
                                               'peptide', 'transitPeptide', 'glycosylation', 'propeptide', 'from',
                                               'to', 'template', 'chain', 'score', 'pdbID', 'pdbSequence', 'fasta'])
        no_info = pd.DataFrame(columns=['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                                        'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence',
                                        'wt_sequence_match', 'whichIsoform', 'datapoint'])
        rest = pd.DataFrame(columns=['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                                     'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence',
                                     'wt_sequence_match', 'whichIsoform', 'datapoint'])

        rest = rest[to_swiss_columns]
        rest = rest.drop_duplicates()

        rest.reset_index(inplace=True)
        rest.drop(['index'], axis=1, inplace=True)
        rest = rest.astype('str')
        to_modbase_size = 0

    print('Modbase matching is completed...\n')
    print('SUMMARY')
    print('-------')
    print('%d data points that failed to match a UniProt Sequence are discarded.' % len(
        not_match_in_uniprot.drop_duplicates(['datapoint'])))
    print('Of the remaining %d:' % uniprot_matched_size)
    print('--%d of %d successfully aligned with PDB structures.' % (
        len(pdb_aligned.drop_duplicates(['datapoint'])), with_pdb_size))
    print('--%d of %d successfully aligned with SwissModels structures.' % (
        len(swiss_match.drop_duplicates(['datapoint'])), to_swiss_size))
    print('--%d of %d successfully aligned with Modbase structures.\n' % (
        len(modbase_match.drop_duplicates(['datapoint'])), to_modbase_size))
    print('--Remaining %d not found to match any models.' % len(rest.drop_duplicates(['datapoint'])))
    print('--A total of %d datapoints will not be evaluated.\n' % (
            len(rest.drop_duplicates(['datapoint'])) + len(not_match_in_uniprot.drop_duplicates(['datapoint']))))

    print('FOR CHECKING : ',
          len(rest.drop_duplicates(['datapoint'])) + len(not_match_in_uniprot.drop_duplicates(['datapoint'])) + len(
              pdb_aligned.drop_duplicates(['datapoint'])) + len(swiss_match.drop_duplicates(['datapoint'])) + len(
              modbase_match.drop_duplicates(['datapoint'])) == data_size)
    no_info = None
    align = None
    not_in_aligned = None
    not_models = None
    modbase_not_match = None



    # Final corrections

    # Now 3D alignment.
    pdb = pdb_aligned.copy()
    swiss = swiss_match.copy()
    modbase = modbase_match.copy()
    pdb_aligned = None
    swiss_match = None
    modbase_match = None

    """
    Fix the columns names for all. 
    WHAT DO WE HAVE NOW?
    - uniprot sequence not found
    - pdb aligned
    - swiss aligned
    - modbase aligned
    - not aligned with anything (rest)
    """

    # Fix the axes and  merge all data.
    # pdb = pd.read_csv('/Users/fatmacankara/Desktop/new_benchmark/pdb_aligned.txt', sep='\t')
    ##pdb = pdb.drop(['pdbInfo', 'resolution'], axis=1)
    # swiss = pd.read_csv('/Users/fatmacankara/Desktop/new_benchmark/swiss_aligned_match.txt', sep=' ')
    # modbase = pd.read_csv('/Users/fatmacankara/Desktop/new_benchmark/modbase_aligned_match.txt', sep='\t')

    pdb.drop(['pdbInfo'], axis=1, inplace=True)
    pdb.rename(columns={'resolution': 'score'}, inplace=True)
    swiss.rename(columns={'qmean_norm': 'score'}, inplace=True)
    modbase.rename(columns={'qmean_norm': 'score'}, inplace=True)
    swiss = swiss[pdb.columns]
    modbase = modbase[pdb.columns]
    pdb['source'] = 'PDB'
    swiss['source'] = 'SWISSMODEL'
    modbase['source'] = 'MODBASE'
    data = pd.concat([swiss, modbase, pdb])

    data.reset_index(inplace=True)
    data.drop(['index'], axis=1, inplace=True)
    data = data.astype('str')
    data_spare = pd.concat([not_match_in_uniprot, rest])
    not_match_in_uniprot = None
    pdb = None
    swiss = None
    modbase = None
    rest = None

    print('Generating FreeSASA files...')
    print('------------------------------------\n')
    data.to_csv(buffer , sep='\t', index=False)
    # Folder to calculated RSA values.

    existing_free_sasa = glob.glob(path_to_output_files + 'freesasa_files/*')
    existing_free_sasa = [i.split('/')[-1].split('.')[0] for i in existing_free_sasa]

    print('Calculation RSA for PDB Structure Files...\n')
    pdb_only = data[data.source == 'PDB']
    for pdbID in pdb_only.pdbID.to_list():
        if pdbID not in existing_free_sasa:
            (run_freesasa(path_to_output_files + 'pdb_structures/' + pdbID.lower() + '.txt',
                          path_to_output_files + 'freesasa_files/' + pdbID.lower() + '.txt', include_hetatms=True,
                          outdir=None, force_rerun=False))


    print('Calculation RSA for SwissModel Files...\n')
    swiss_only = data[data.source == 'SWISSMODEL']
    swiss_dp = []
    for i in swiss_only.index:
        swiss_dp.append(swiss_only.at[i, 'uniprotID'] + '_' + swiss_only.at[i, 'pdbID'].lower() + '_' + str(
            round(float(swiss_only.at[i, 'score']), 2)))
    for pdbID in swiss_dp:
        if pdbID not in existing_free_sasa:
            (run_freesasa(path_to_output_files + 'swissmodel_structures/' + pdbID + '.txt',
                          path_to_output_files + 'freesasa_files/' + pdbID + '.txt', include_hetatms=True,
                          outdir=None, force_rerun=False))

    print('Calculation RSA for Modbase Model Files...\n')
    modbase_only = data[data.source == 'MODBASE']
    for pdbID in modbase_only.pdbID.to_list():
        if pdbID not in existing_free_sasa:
            (run_freesasa(path_to_output_files + 'modbase_structures_individual/' + pdbID.lower() + '.txt',
                          path_to_output_files + 'freesasa_files/' + pdbID.lower() + '.txt', include_hetatms=True,
                          outdir=None, force_rerun=False))

    # This annotation list is different than the prev one, keep it.
    annotation_list = ['disulfide', 'intMet', 'intramembrane', 'naturalVariant', 'dnaBinding', 'activeSite',
                       'nucleotideBinding', 'lipidation', 'site', 'transmembrane', 'crosslink', 'mutagenesis', 'strand',
                       'turn', 'helix', 'metalBinding', 'repeat', 'caBinding', 'topologicalDomain', 'bindingSite',
                       'region',
                       'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil', 'peptide',
                       'transitPeptide', 'glycosylation', 'propeptide', 'domainStartonPDB', 'domainEndonPDB']

    folder_path = path_to_output_files + 'freesasa_files/'

    aligner = Align.PairwiseAligner()
    print('Proceeding to 3D distance calculation...\n')

    data.domainEndonPDB = data.domainEndonPDB.astype(str)
    data.domainStartonPDB = data.domainStartonPDB.astype(str)

    existing_free_sasa = None
    swiss_dp = None
    pdb_only = None
    swiss_only = None
    modbase_only = None

    for i in data.index:
        print(i,'of', len(data), data.at[i,'uniprotID'])
        if data.at[i, 'source'] == 'PDB':
            pdb_path = path_to_output_files + 'pdb_structures/' + data.at[i, 'pdbID'] + '.txt'
        elif data.at[i, 'source'] == 'MODBASE':
            pdb_path = path_to_output_files + 'modbase_structures_individual/' + data.at[i, 'pdbID'] + '.txt'
        elif data.at[i, 'source'] == 'SWISSMODEL':
            pdb_path = path_to_output_files + 'swissmodel_structures/' + data.at[i, 'uniprotID'] + '_' + data.at[
                i, 'pdbID'] + '_' + str(data.at[i, 'score']) + '.txt'
        data['uniprotSequence'] = data['uniprotSequence'].str.replace('U', 'C')
        data['pdbSequence'] = data['pdbSequence'].str.replace('U', 'C')
        alignments = get_alignments_3D(data, pdb_path)
        mutPos = data.at[i, 'mutationPositionOnPDB']
        try:
            coordMut = get_coords(mutPos, alignments)[0]
        except:
            ValueError
            coordMut = 'nan'
        try:
            sasa_pos = get_coords(mutPos, alignments)[2]
            data.at[i, 'sasa'] = sasa(data, sasa_pos)
        except:
            ValueError
            data.at[i, 'sasa'] = 'nan'  # mutation position is nan
        for annot in annotation_list:
            annotx = []
            try:
                positions_of_annotations = data.at[i, annot].split(',')
                for pos in positions_of_annotations:
                    pos = pos.strip().strip('\'').strip('[\'').strip('\']')
                    try:
                        if '-' not in pos:
                            pos = int(float(pos))
                            coordAnnot = get_coords(pos, alignments)[0]
                            try:
                                annotx.append(find_distance(coordMut, coordAnnot))
                            except:
                                ValueError  # Herhangi bir distance ekleme

                        else:
                            for r in range(int(pos.split('-')[0]), int(pos.split('-')[1]) + 1):
                                coordAnnot = get_coords(r, alignments)[0]
                                annotx.append(find_distance(coordMut, coordAnnot))
                    except:
                        ValueError
                try:
                    data.at[i, annot] = min([float(i) for i in annotx])
                except:
                    ValueError
                    data.at[i, annot] = 'nan'

            except:
                ValueError

        if (str(data.at[i, 'domainStartonPDB']) == 'NaN' or str(data.at[i, 'domainStartonPDB']) == 'nan') and (
                str(data.at[i, 'domainEndonPDB']) != 'NaN' and str(data.at[i, 'domainEndonPDB']) != 'nan'):
            data.at[i, 'domainStartonPDB'] = 100000
        elif (str(data.at[i, 'domainEndonPDB']) == 'NaN' or str(data.at[i, 'domainEndonPDB']) == 'nan') and (
                str(data.at[i, 'domainStartonPDB']) != 'NaN' and str(data.at[i, 'domainStartonPDB']) != 'nan'):
            data.at[i, 'domainEndonPDB'] = 100000
        elif (str(data.at[i, 'domainStartonPDB']) == 'NaN' and str(data.at[i, 'domainEndonPDB']) == 'nan'):
            data.at[i, 'domaindistance3D'] = 'nan'

        data.at[i, 'domaindistance3D'] = min(float(data.at[i, 'domainStartonPDB']),
                                             float(data.at[i, 'domainEndonPDB']))
        data.at[i, 'domaindistance3D'] = min(float(data.at[i, 'domainStartonPDB']),
                                             float(data.at[i, 'domainEndonPDB']))
        # (data[data.datapoint == data.at[i, 'datapoint']]).astype(str).to_csv(path_to_output_files + 'after_3D_alignment_append.txt',
        #                                                                     mode='a', index=False, sep='\t', header=False)

    data = data.astype(str)
    data.replace({'NaN': 'nan'}, inplace=True)

    # data.to_csv(path_to_output_files + 'after_3D_alignment.txt', sep='\t', index=False)

    # Now unify all 3 separate data. We have with_pdb. The ones that have pdb structyres, swiss, modbase, the ones didnt match with ant and the ones didnt have wt seq match.

    # Get interface positions from ECLAIR. Download HQ human
    print()
    print('Assigning surface regions...')
    print('------------------------------------\n')

    print('Extracting interface residues...\n')
    data_interface = pd.read_csv(path_to_interfaces, sep='\t')

    positions = get_interface_positions(data_interface, 'P1', 'P2')

    interface_dataframe = pd.DataFrame()

    for key, val in positions.items():
        k = pd.Series((key, str(list(set(val)))))
        interface_dataframe = interface_dataframe.append(k, ignore_index=True)
    interface_dataframe.columns = ['uniprotID', 'positions']

    if len(data) == 0:
        data = pd.DataFrame(columns=['uniprotID', 'wt', 'mut', 'pos', 'composition', 'polarity', 'volume',
                                     'domain', 'domStart', 'domEnd', 'distance', 'uniprotSequence',
                                     'pdbSequence', 'wt_sequence_match', 'whichIsoform', 'pdbID', 'score',
                                     'chain', 'datapoint', 'disulfide', 'intMet', 'intramembrane',
                                     'naturalVariant', 'dnaBinding', 'activeSite', 'nucleotideBinding',
                                     'lipidation', 'site', 'transmembrane', 'crosslink', 'mutagenesis',
                                     'strand', 'helix', 'turn', 'metalBinding', 'repeat',
                                     'topologicalDomain', 'caBinding', 'bindingSite', 'region',
                                     'signalPeptide', 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil',
                                     'peptide', 'transitPeptide', 'glycosylation', 'propeptide',
                                     'disulfideBinary', 'intMetBinary', 'intramembraneBinary',
                                     'naturalVariantBinary', 'dnaBindingBinary', 'activeSiteBinary',
                                     'nucleotideBindingBinary', 'lipidationBinary', 'siteBinary',
                                     'transmembraneBinary', 'crosslinkBinary', 'mutagenesisBinary',
                                     'strandBinary', 'helixBinary', 'turnBinary', 'metalBindingBinary',
                                     'repeatBinary', 'topologicalDomainBinary', 'caBindingBinary',
                                     'bindingSiteBinary', 'regionBinary', 'signalPeptideBinary',
                                     'modifiedResidueBinary', 'zincFingerBinary', 'motifBinary',
                                     'coiledCoilBinary', 'peptideBinary', 'transitPeptideBinary',
                                     'glycosylationBinary', 'propeptideBinary', 'pdb_alignStatus',
                                     'mutationPositionOnPDB', 'domainStartonPDB', 'domainEndonPDB',
                                     'source', 'sasa', 'domaindistance3D', 'threeState_trsh4_HQ', 'domain_fisher'])
    else:
        data.sasa = data.sasa.astype('str')

    for i in data.index:
        if '*' in data.at[i, 'sasa']:
            data.at[i, 'sasa'] = data.at[i, 'sasa'].split('*')[0]

    data.sasa = data.sasa.replace({'N/A': 'nan'})
    data.sasa = data.sasa.replace({'None': 'nan'})
    data.replace({'   N/A': 'nan'}, inplace=True)
    data.replace({'None': 'nan'}, inplace=True)
    # print('Without calculated rsa (nan):', len(data[data.sasa == 'nan']))
    # print('With calculated rsa (nan):', len(data[data.sasa != 'nan']))
    data.sasa = data.sasa.astype(float)
    data = data.astype(str)
    for i in data.index:
        if float(data.at[i, 'sasa']) < 5:
            data.at[i, 'trsh4'] = 'core'
        elif float(data.at[i, 'sasa']) >= 5:
            data.at[i, 'trsh4'] = 'surface'
        elif data.at[i, 'sasa'] == 'nan':
            data.at[i, 'trsh4'] = 'nan'

    data = data.merge(interface_dataframe, on='uniprotID', how='left')
    data.positions = data.positions.astype('str')
    for i in data.index:
        if (str(data.at[i, 'pos']) in data.at[i, 'positions']) and data.at[i, 'trsh4'] == 'surface':
            print((str(data.at[i, 'pos']) in data.at[i, 'positions']))
            data.at[i, 'threeState_trsh4_HQ'] = 'interface'
        elif (str(data.at[i, 'pos']) not in data.at[i, 'positions']) and data.at[i, 'trsh4'] == 'surface':
            data.at[i, 'threeState_trsh4_HQ'] = 'surface'
        elif (str(data.at[i, 'pos']) not in data.at[i, 'positions']) and data.at[i, 'trsh4'] == 'core':
            data.at[i, 'threeState_trsh4_HQ'] = 'core'
        elif (str(data.at[i, 'pos']) in data.at[i, 'positions']) and data.at[i, 'trsh4'] == 'core':
            data.at[i, 'threeState_trsh4_HQ'] = 'conflict'
        elif data.at[i, 'trsh4'] == 'nan':
            data.at[i, 'threeState_trsh4_HQ'] = 'nan'

    data.drop(['positions'], axis=1, inplace=True)

    """
    print('Counts for HQ, treshold 4:')
    print('HQ_4_surface: ', len(data[data.threeState_trsh4_HQ == 'surface']))
    print('HQ_4_core: ', len(data[data.threeState_trsh4_HQ == 'core']))
    print('HQ_4_interface: ', len(data[data.threeState_trsh4_HQ == 'interface']))
    print('HQ_4_conflict: ', len(data[data.threeState_trsh4_HQ == 'conflict']))
    print(len(data[data.threeState_trsh4_HQ == 'surface']) + len(data[data.threeState_trsh4_HQ == 'interface']),
          'should be equal to', len(finalData[finalData.trsh4 == 'surface']))
    print(len(data[data.threeState_trsh4_HQ == 'core']) + len(data[data.threeState_trsh4_HQ == 'conflict']), 'should be equal to',
          len(data[data.trsh4 == 'core']))

    HQ_4_surface:  202
    HQ_4_core:  256
    HQ_4_interface:  90
    HQ_4_conflict:  35
    """

    # OPTIONAL
    # DOMAIN SELECTION

    # Next step: Delete all other domains with 'Domain X.' R is capable of handling 53 categories. We will keep 52 most
    # significant domains and 53th category will be Domain X.
    #
    fisherResult = pd.read_csv(fisher_path, sep='\t')

    significant_domains = fisherResult.domain.to_list()
    for i in data.index:
        if data.at[i, 'domain'] in significant_domains:
            data.at[i, 'domain_fisher'] = data.at[i, 'domain']
        else:
            data.at[i, 'domain_fisher'] = 'domainX'

    # Change the numbering for binary annotations and create 3 classes:
    # nan--> 0, 0 -->1 and 1 -->2

    print('Final adjustments are being done...\n')
    binaryCols = ['disulfideBinary', 'intMetBinary', 'intramembraneBinary', 'naturalVariantBinary', 'dnaBindingBinary',
                  'activeSiteBinary', 'nucleotideBindingBinary', 'lipidationBinary', 'siteBinary',
                  'transmembraneBinary', 'crosslinkBinary', 'mutagenesisBinary',
                  'strandBinary', 'helixBinary', 'turnBinary', 'metalBindingBinary',
                  'repeatBinary', 'caBindingBinary', 'topologicalDomainBinary',
                  'bindingSiteBinary', 'regionBinary', 'signalPeptideBinary',
                  'modifiedResidueBinary', 'zincFingerBinary', 'motifBinary',
                  'coiledCoilBinary', 'peptideBinary', 'transitPeptideBinary',
                  'glycosylationBinary', 'propeptideBinary']
    data = data.astype(str)
    data.replace({'NaN': 'nan'}, inplace=True)
    for i in data.index:
        for j in binaryCols:
            data[j] = data[j].astype('str')
            if (data.at[i, j] == '0') or (data.at[i, j] == '0.0'):
                data.at[i, j] = '1'
            elif data.at[i, j] == 'nan':
                data.at[i, j] = '0'
            elif (data.at[i, j] == '1') or (data.at[i, j] == '1.0'):
                data.at[i, j] = '2'

    annotCols = ['disulfide', 'intMet', 'intramembrane',
                 'naturalVariant', 'dnaBinding', 'activeSite', 'nucleotideBinding',
                 'lipidation', 'site', 'transmembrane', 'crosslink', 'mutagenesis',
                 'strand', 'helix', 'turn', 'metalBinding', 'repeat', 'caBinding',
                 'topologicalDomain', 'bindingSite', 'region', 'signalPeptide',
                 'modifiedResidue', 'zincFinger', 'motif', 'coiledCoil', 'peptide',
                 'transitPeptide', 'glycosylation', 'propeptide']
    # Bu kod binary 2 olup 3ddistanceı 0dan farklı olanların oranını bulmak için.
    """
    print('annotationName', 'class2_distanceNot0', 'class2_distance0', 'ratio in total')
    print('--------------' * 3)
    for annot in annotCols:
        countWrong = 0
        countRight = 0
        binaryName = str(annot) + 'Binary'
        for i in data.index:
            if data.at[i, binaryName] == '2' and (data.at[i, annot] != '0.0' and data.at[i, annot] != '0') and (data.at[i, annot] != 'nan'):
                countWrong +=1
            elif data.at[i, binaryName] == '2' and (data.at[i, annot] == '0.0' or data.at[i, annot] == '0'):
                countRight += 1
        print(annot, countWrong, countRight)
    """
    for i in data.index:
        for annot in annotCols:
            binaryName = str(annot) + 'Binary'
            if data.at[i, binaryName] == '2':
                data.at[i, annot] = '0.0'
    data.replace({'100000': 'nan'}, inplace=True)
    # data.to_csv(path_to_output_files + 'benchmark_completed4.txt', sep='\t', index=False)
    data = add_physicochemical(data)
    final_cols_to_keep = ['datapoint', 'composition', 'polarity', 'volume', 'granthamScore','domain', 'domain_fisher',
                          'domaindistance3D', 'disulfide', 'intMet', 'intramembrane', 'naturalVariant',
                          'dnaBinding', 'activeSite', 'nucleotideBinding', 'lipidation', 'site',
                          'transmembrane', 'crosslink', 'mutagenesis', 'strand', 'helix', 'turn',
                          'metalBinding', 'repeat', 'topologicalDomain', 'caBinding',
                          'bindingSite', 'region', 'signalPeptide', 'modifiedResidue',
                          'zincFinger', 'motif', 'coiledCoil', 'peptide', 'transitPeptide',
                          'glycosylation', 'propeptide', 'disulfideBinary', 'intMetBinary',
                          'intramembraneBinary', 'naturalVariantBinary', 'dnaBindingBinary',
                          'activeSiteBinary', 'nucleotideBindingBinary', 'lipidationBinary',
                          'siteBinary', 'transmembraneBinary', 'crosslinkBinary',
                          'mutagenesisBinary', 'strandBinary', 'helixBinary', 'turnBinary',
                          'metalBindingBinary', 'repeatBinary', 'topologicalDomainBinary',
                          'caBindingBinary', 'bindingSiteBinary', 'regionBinary',
                          'signalPeptideBinary', 'modifiedResidueBinary', 'zincFingerBinary',
                          'motifBinary', 'coiledCoilBinary', 'peptideBinary',
                          'transitPeptideBinary', 'glycosylationBinary', 'propeptideBinary', 'sasa',
                          'threeState_trsh4_HQ']

    ready = data[final_cols_to_keep]
    ready.to_csv(outpath, sep='\t', index=False)

    print('Feature vector successfully created...')

end = timer()
hours, rem = divmod(end - start, 3600)
minutes, seconds = divmod(rem, 60)
print("Time passed: {:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds))
sys.stdout.close()
