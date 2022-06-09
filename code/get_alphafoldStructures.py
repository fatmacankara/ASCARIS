import tarfile, glob, os
from biopandas.pdb import PandasPdb
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='ASCARIS')

parser.add_argument('-file_name', '--file_name',
                    help='Enter the file tar file name to untar',
                    default=1)

args = parser.parse_args()

alphafold = args.file_name

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
# Unzip AlphaFold structures

def create_file():
    os.makedirs('input_files/alphafold_structures/', exist_ok=True)
    for f in glob.glob(f'input_files/{alphafold}'):
        with tarfile.open(f) as tar:
            tar.extractall(f'input_files/alphafold_structures/')

    # Create summary file
    alphafold_summary_file = open('input_files/alphafold_summary.txt', 'w')
    alphafold_summary_file.write('uniprotID\tchain\tsequence\tmodel_num')
    alphafold_summary_file.write('\n')
    for f in glob.glob('input_files/alphafold_structures/*pdb*'):
        str1 = PandasPdb().read_pdb(f)
        str1 = str1.df['ATOM']
        str1 = str1[['alt_loc', 'residue_name', 'residue_number', 'atom_name', 'insertion', 'chain_id']]
        str1 = str1[str1.atom_name == 'CA']
        str1['residue_name'] = str1['residue_name'].apply(lambda x: threeToOne(x))
        str1['alt_loc'] = str1['alt_loc'].replace({'': np.NaN})
        str1 = str1.drop_duplicates(['residue_name', 'residue_number'])
        structure_residues_pdb = ''.join(str1.residue_name.to_list())
        model_no = f.split('-')[2].strip()[1:]
        up_name = f.split('-')[1].strip()
        chain_id = list(set(str1.chain_id.to_list()))[0]
        alphafold_summary_file.write(up_name)
        alphafold_summary_file.write('\t')
        alphafold_summary_file.write(chain_id)
        alphafold_summary_file.write('\t')
        alphafold_summary_file.write(structure_residues_pdb)
        alphafold_summary_file.write('\t')
        alphafold_summary_file.write(model_no)
        alphafold_summary_file.write('\n')


if __name__ == '__main__':
    create_file()