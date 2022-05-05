import glob
import ssbio.utils
import subprocess
import ssbio
import os.path as op
from add_3Dalignment import *
import os

import gzip
import shutil



def run_freesasa(infile, outfile, include_hetatms=True, outdir=None, force_rerun=False, file_type = 'gzip'):
    if not outdir:
        outdir = ''
    outfile = op.join(outdir, outfile)
    if file_type == 'pdb':
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
    elif file_type == 'gzip':
        with gzip.open(infile, 'rb') as f_in:
            with open('file_temp.pdb', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        infile = 'file_temp.pdb'

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
    os.remove(infile)
    return outfile

def calculate_freesasa(ID, model_num, existing_free_sasa, path_to_input,path_to_output_files, file_type = 'gzip'):
    print('Calculating surface area...\n')
    if file_type == 'gzip':
        if ID not in existing_free_sasa:
            fullID = f'AF-{ID}-F{model_num}-model_v1.pdb.gz'
            run_freesasa(f'{path_to_input}{fullID}',  # BURDAKİ INPUT PATHİ DEĞİŞMELİ, NORMALDE EXAMPLE INUT BBNUN İÇİNDE DEĞİL. PATH TO OUT YERİNE BİR ÜST FOLDE KULLANILANİLİR.
                         f'{path_to_output_files}/freesasa_files/{fullID}.txt', include_hetatms=True,
                         outdir=None, force_rerun=False)
    elif file_type == 'pdb':
        if ID not in existing_free_sasa:
            fullID = f'AF-{ID}-F{model_num}-model_v1.pdb'
            run_freesasa(f'{path_to_input}{fullID}',  # BURDAKİ INPUT PATHİ DEĞİŞMELİ, NORMALDE EXAMPLE INUT BBNUN İÇİNDE DEĞİL. PATH TO OUT YERİNE BİR ÜST FOLDE KULLANILANİLİR.
                         f'{path_to_output_files}/freesasa_files/{fullID}.txt', include_hetatms=True,
                         outdir=None, force_rerun=False)

def sasa(source, pdbID, uniprotID, sasa_pos, wt, mode, path_to_output_files,file_type = 'gzip'):
    if mode == 1:
        sasa = 'nan'
        for filename in glob.glob(f'{path_to_output_files}/freesasa_files/*'):
            if source == 'PDB':
                fname = filename.split('.')[0].split('/')[-1].upper()
            elif source == 'MODBASE':
                fname = filename.split('.')[0].split('/')[-1]
            elif source == 'SWISSSMODEL':
                fname = filename.split('_')[2]

            if pdbID == fname:
                files = open(filename, 'r')
                file = files.readlines()
                for k in file:
                    if k.strip()[10:13] == sasa_pos:
                        residue = str(k[4:7].strip())
                        if wt == threeToOne(residue):
                            sasa = str(k[22:28]).strip('\n')
                            return (sasa)
                        elif wt != threeToOne(residue):
                            sasa = str(k[22:28]).strip('\n') + '*'
                            return (sasa)
                        else:
                            return 'nan'  #######

    if mode == 2:
        if sasa_pos != np.NaN:
            sasa = 'nan'
            if file_type == 'pdb':
                for filename in glob.glob(f'{path_to_output_files}/freesasa_files/*'):
                    fname = list(filter(None, filename.split('.'))).split('/')[-1].upper()
                    if uniprotID == fname:
                        files = open(filename, 'r')
                        file = files.readlines()
                        for k in file:
                            if k.strip()[10:13] == sasa_pos:
                                residue = str(k[4:7].strip())
                                if wt == threeToOne(residue):
                                    sasa = str(k[22:28]).strip('\n')
                                elif wt != threeToOne(residue):
                                    sasa = str(k[22:28]).strip('\n') + '*'

                return sasa  #######
            elif file_type == 'gzip':
                for filename in glob.glob(f'{path_to_output_files}/freesasa_files/*'):
                    fname = list(filter(None, filename.split('.')))[0].split('/')[-1].split('-')[1].upper()

                    if uniprotID == fname:
                        files = open(filename, 'r')
                        file = files.readlines()
                        for k in file:
                            if str(k.strip()[10:13]) == str(sasa_pos):
                                residue = str(k[4:7].strip())
                                if wt == threeToOne(residue):
                                    sasa = str(k[22:28]).strip('\n')
                                elif wt != threeToOne(residue):
                                    sasa = str(k[22:28]).strip('\n') + '*'
                                else:
                                    sasa = 'nan'

                return sasa
        else:
            sasa = 'nan'
            return sasa