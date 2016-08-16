'''
------------------------------------------------------------------------------
AUTHOR: Raphael Peer, raphael1peer@gmail.com

PURPOSE:
From a directory containing PDB-files, extract all sequences (from the actual
structures, not the header) and write them to a fasta file.
Apply the script to the structures which will be used for residue contact
calculation. The extracted sequences can ee used to create a reference
alignment with a sequence-alignment tool, such as MUSCLE or ClustalW.


INPUT:
1) Directory with PDB-files:
Ony PDB-files with exactly one chain are accepted. The reason for this is that
the script was created to extract the seqeunces from the PDB-files used for
residue conctact calculation. These PDB-files are pre-processed to have only
one chain per file which contains the domain of interest (see README or
docstring of 'process_pdb.py'). It is recommended to first execute the script 
'process_pdb.py'. This script prepares 'raw' PDB-files for residue contact
calculation and puts them in the newly created directory 'processed_pdb_files'.
Provide the path to this directory as first argument.
2) Name of the output file (optional).

NOTE:
There are two ways of creating a reference alignment.
1) Use structure alignment (e.g. with MUSTANG). Obviosly, the exact same
structures used for residue contact calculation have to be used for the
structural alignment (i.e. pre-processed with 'process_pdb.py')
2) Use sequence alignemnt (e.g. with MUSCLE or ClustalW). This script can be
used to extract the sequences from the PDB-files of interest. Then, the
extracted sequences can be aligned to create the reference alignment.
------------------------------------------------------------------------------
'''

import os
import sys
import argparse
from Bio import PDB


parser = argparse.ArgumentParser()
parser.add_argument('processed_pdb_dir', help='Directory with PDB-structures. '
                    'Only structure with exactly one chain will be accepted. '
                    'It is recommended to run the script "process_pdb.py" and '
                    'use the output ("processed_pdb_files") as input here.')
parser.add_argument('outfile_name', nargs='?',type=str,
                    const='PDB_sequences.fa', default='PDB_sequences.fa',
                    help='Name of the output file (optional)')
try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(1)


def write_pdb_seq_to_file(pdb_file):
    '''IN: (path to) PDB-file with only one chain
    OUT: sequence of PDB-file (from actual structure, not header)'''
    struct = PDB.PDBParser().get_structure('current', pdb_file)
    assert len(list(struct.get_chains())) == 1, \
        'WARINING: There are more than one chains in structure %s. \
        \n It will be excluded from analysis.' % struct.get_id()
    seq = ''
    for pp in PDB.PPBuilder().build_peptides(struct):
        seq += pp.get_sequence()
    return seq


with open(args.outfile_name, 'w') as outfile:
    for filename in os.listdir(args.processed_pdb_dir):
        sequence = write_pdb_seq_to_file(os.path.join(args.processed_pdb_dir,
                                                      filename))
        outfile.write('>%s\n%s\n' % (filename, sequence))
