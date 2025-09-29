#!/usr/bin/env python3

'''
Script to test the multi molfile_to_params function
on FADE hit ligand candidates using Explorer:
'''

import subprocess
import argparse
import glob
import sys
import os

script_path = '/projects/SimBioSys/share/software/rosetta/rosetta.binary.linux.release-296/main/source/scripts/python/public/molfile_to_params.py'

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('pattern', type=str, help='Glob pattern to use for finding conformer files in sdf format.')
    
    args = parser.parse_args()

    conf_files = sorted(glob.glob(args.pattern))

    for i, file in enumerate(conf_files):
        lig_pdb_tag = f'NM{i+1}'
        prefix = os.path.basename(os.path.splitext(file)[0])
        cmd = [
            sys.executable, script_path,
            '-n', lig_pdb_tag,
            '-p', prefix,
            '--chain=X',
            '--conformers-in-one-file', file
        ]
        subprocess.run(cmd, check=True)

if __name__ == '__main__':
    main()