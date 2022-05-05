import pdb_featureVector
import alphafold_featureVector


input_set = input('Enter Query DataPoint\n'
                  'Option 1: Comma-separated list of idenfiers (UniProt ID-wt residue-position-mutated residue (e.g. Q9Y4W6-N-432-T or Q9Y4W6-N-432-T, Q9Y4W6-N-432-T)) \n'
                  'Option 2: Enter comma-separated file path \n')

mode = int(input('Enter 1 for PDB-Modbase-SwissModel structures, 2 for AlphaFold structures: \n'))

def run_featureVector():
    print('*****************************************')
    print('Feature vector generation is in progress. \nPlease check log file for updates..')
    print('*****************************************')

    if mode == 1:
        pdb_featureVector.pdb(input_set, mode)
    elif mode == 2:
        alphafold_featureVector.alphafold(input_set, mode)

if __name__ == '__main__':
    run_featureVector()


