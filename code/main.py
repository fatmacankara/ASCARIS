import pdb_featureVector
import alphafold_featureVector
import argparse

parser = argparse.ArgumentParser(description='ASCARIS')

parser.add_argument('-o', '--input_option',
                    help='Selection of input structure data.\n 1: PDB Structures (default), 2: AlphaFold Structures',
                    default=1)
parser.add_argument('-i', '--input_datapoint',
                    help='Input file or query datapoint\n Option 1: Comma-separated list of idenfiers (UniProt ID-wt residue-position-mutated residue (e.g. Q9Y4W6-N-432-T or Q9Y4W6-N-432-T, Q9Y4W6-N-432-T)) \n Option 2: Enter comma-separated file path')

args = parser.parse_args()

mode = args.input_option
input_set = args.input_datapoint
print(input_set)
"""
input_set = input('Enter Query DataPoint\n'
                  'Option 1: Comma-separated list of idenfiers (UniProt ID-wt residue-position-mutated residue (e.g. Q9Y4W6-N-432-T or Q9Y4W6-N-432-T, Q9Y4W6-N-432-T)) \n'
                  'Option 2: Enter comma-separated file path \n')

mode = int(input('Enter 1 for PDB-Modbase-SwissModel structures, 2 for AlphaFold structures: \n'))
"""
def run_featureVector(mode, input_set):

    print('*****************************************')
    print('Feature vector generation is in progress. \nPlease check log file for updates..')
    print('*****************************************')
    mode = int(mode)
    if mode == 1:
        pdb_featureVector.pdb(input_set, mode)
    elif mode == 2:
        alphafold_featureVector.alphafold(input_set, mode)

if __name__ == '__main__':
    run_featureVector(mode, input_set)


