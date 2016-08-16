'''
------------------------------------------------------------------------------
AUTHOR: Raphael Peer, raphael1peer@gmail.com

PURPOSE:
The script prepares the input PDB-files for residue contact calculation.

OUTPUT:
Processed PDB-files in new directory (inside the directory 'results')

PDB-FILE PREPARATION:
1) Selection of the relevant chain which contains the Pfam-domain of interest.
If the Pfam-domain is not present in the structure, no chain is selected.
If the Pfam-domain is present on multiple chains, the first chain
(alphabetically) is selected.
2) Removal of all lines which do not start with 'ATOM', 'TER' or 'END'. This
includes heteroatoms (HETATM) or anisotropic temperature factors (ANISOU)
3) Removal of hydrogen atoms (if present)

NOTE:
The fully automated large scale analysis of protein structures requires
standardized input data. For instance, most structures in the PDB do not include hydrogen
atoms (as they were solved with x-ray crystallography). For this reason,
hydrogen atoms are removed if present to allow comparison with other
structures.
------------------------------------------------------------------------------
'''

import argparse
import sys
import os
import pandas as pd
from Bio import PDB

parser = argparse.ArgumentParser()
parser.add_argument('raw_pdb_dir', help='directory with raw (unprocessed) '
                    'PDB-structures')
parser.add_argument('pfam_domain', help='Pfam domain of interest. Important '
                    'for the program to know which chains to analyse in '
                    'case of complex structures.')
parser.add_argument('sifts_chain_pfam', help='SIFTS-file "pdb_chain_pfam.csv" for '
                    'finding chains withing the PDB-structures which contain '
                    'the Pfam-domain of interest')
parser.add_argument('processed_pdb_dir', help='Name of the output directory for the '
                    'processed PDB-files.')
try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(1)
chain_pfam = pd.read_csv(args.sifts_chain_pfam, comment='#')


# CHECK IF OUPUT DIRECTORY EXISTS
if not os.path.exists(args.processed_pdb_dir):
    os.makedirs(args.processed_pdb_dir)
elif len(os.listdir(args.processed_pdb_dir)) > 0:
    print 'Directory %s is not empty.' % args.processed_pdb_dir
    user_input =  raw_input('Do you want to delete its content?\n(y for yes) ')
    if user_input == 'y':
        os.system('rm -f %s/*' % args.processed_pdb_dir)
    else:
        print('Please make sure the directory is empty or provide a different '
            'name for the output directory. Then restart script.')
        sys.exit(1)


def select_chain(filename, domain_of_int):
    '''selects only one chain per structure; chain has to have domain of
    interest (and this domain only'''
    pdb_id = filename.split('.')[0]
    sub = chain_pfam[(chain_pfam.PDB == pdb_id) & (chain_pfam.PFAM_ID == domain_of_int)]
    chains_of_interest = list(sub.CHAIN)
#    print '\n%s' % pdb_id
#    print 'all chains: %s' % list(chain_pfam[(chain_pfam.PDB == pdb_id)].CHAIN)
#    print 'chains of interest: %s' % chains_of_interest
    if len(chains_of_interest) >= 1:
        return chains_of_interest[0]   # return first chain with domain of interest
# 4mit and 1u8z are not included in the chain_pfam file.
# Hence no chains are selected

def extract_chain(in_dir, filename, chain_of_int, out_dir):
    '''saves chain of interest as separate structure'''
    pdb_id = filename.split('.')[0]
    struct = PDB.PDBParser().get_structure(pdb_id, os.path.join(in_dir, filename))
    io = PDB.PDBIO()
    io.set_structure(struct[0][chain_of_int])
    io.save(os.path.join(out_dir, filename))


def clean(directory, filename):
    '''greps only ATOM, TER or END-lines from pdb, overwrites original file'''
    if not os.path.exists(os.path.join(directory, filename)):
        print 'file %s not found. Script continues.'
        return
    os.system('mv %s tmp.pdb' % os.path.join(directory, filename))  # rename original file
    os.system('egrep "^ATOM|^TER|^END" tmp.pdb > %s' % os.path.join(directory, filename))
    os.system('rm tmp.pdb')


def remove_hydrogens(directory, filename):
    '''Removes all hydrogen atoms from a PDB-file, overwrites original file
    Note: requires "cleaned" input PDB-file (as it ignores HETATM)'''
    if not os.path.exists(os.path.join(directory, filename)):
        print 'file %s not found. Script continues.'
        return
    os.system('mv %s tmp.pdb' % os.path.join(directory, filename))  # rename original file
    os.system('cat tmp.pdb | sed "s/[ \t]*$//" | grep -v "^ATOM.*H$" > %s' % os.path.join(directory, filename))
    os.system('rm tmp.pdb')


selected_chains = pd.DataFrame(columns = ['pdb_id', 'chain'])  # just for info
print('\nPDB-ID\tchain selected for analysis')
for pdb_file in os.listdir(args.raw_pdb_dir):
    chain = select_chain(pdb_file, args.pfam_domain)
    print "%s\t%s" % (pdb_file, chain)
    if chain != None:  # if chain with Pfam-domain of interest found
        selected_chains = selected_chains.append({'pdb_id': pdb_file.split('.')[0],
                                                  'chain': chain},
                                                  ignore_index = True)
        extract_chain(args.raw_pdb_dir, pdb_file, chain, args.processed_pdb_dir)
        clean(args.processed_pdb_dir, pdb_file)
        remove_hydrogens(args.processed_pdb_dir, pdb_file)

selected_chains.to_csv('results/selected_chains_info.csv', index=False)

print '\nOut of %s raw PDB-files, for %s PDB-files, a chain with the Pfam \
domain of interst could be extracted and written to a new file in the \
directory %s' % (len(os.listdir(args.raw_pdb_dir)), len(selected_chains),
                 args.processed_pdb_dir)
   
