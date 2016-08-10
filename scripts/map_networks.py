'''
------------------------------------------------------------------------------
AUTHOR: Raphael Peer, raphael1peer@gmail.com

PURPOSE:
1) For each residue of each structure, identify the position in the reference
alignment.
2) For each residue of each structure, identify the equivalent residue in a
chosen reference structure.

OUTPUT:
csv-file: PDB-ID, residue number (always starting at 1), PDB-residue number
(as given in the PDB-file), alignment position number, amino acid, PDB-residue
number of the reference structure
Example:
pdb_id,resnum,pdb,alignment_pos,aa,ref_pdb
1gwn,2.0,23.0,16.0,LYS,19.0
The second residue in 1gwn.pdb (a lysine) has the number 23 in the PDB-file.
In the reference alignment, this residue is at position 16.
The equivalent residue in the reference structure (in this case 1g16) has
the number 19 in 1g16.pdb. 'Structurally equivalent residue' means that, when
the structures are superimposed in structural alignment, the two residues are
at the same place. In the reference alignment, created from the structural
alignment (in this case with the software MUSTANG), those 'equivalent residues'
are in the same column.

NOTE1:
The reference structure can be any structure from the dataset and is provided
as an argument to this script. The feature allows to interpret the results
using a structure of particular interest. However, if the feature is not
needed, any random structure from the dataset can be defined as reference-
structure and the column 'ref_pdb' in the results can simply be ignored.

NOTE2:
The datatype of the residue numbers in the output csv-file are provided as
floats rather than integers as one would intiuitively expect. This is due to
the fact, that there is no 'NA' for integers in pandas. As 'NA'-values may
(if there is no structurally equivalent residue in the reference structure),
residue numbers cannot be converted to integers.
------------------------------------------------------------------------------
'''

import argparse
import sys
import os
import pandas as pd
from Bio import SeqIO
from Bio import PDB

parser = argparse.ArgumentParser()
parser.add_argument('processed_pdb_dir', help='Directory with processed '
                    'PDB-structures for calculation of residue contact '
                    'network')
parser.add_argument('reference_alignment', help='Alignment of the sequences '
                    'of the structures for which residue contact networks '
                    'were created. See the documentation for information '
                    'about the requirements of such an alignment')
parser.add_argument('reference_structure', help='For all residues, the '
                    'equivalent residues (PDB-numbering) of the reference '
                    'structure will be provided. Just provide the PDB-ID of '
                    'your favourite structure of the dataset')
try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(1)


def map_to_alignment(alignment, sequence_id):
    '''IN: multiple sequence alingment (fasta), identifier of one sequence
    OUT: referenc list: residue number as index, alignment position as value'''
    id2seq = SeqIO.to_dict(SeqIO.parse(alignment, 'fasta'))
    alignment_positions = list(id2seq[sequence_id])
    ref_list = []
    for i in range(len(alignment_positions)):
        if alignment_positions[i] != '-':
            ref_list.append(i + 1)
    return(ref_list)


def map_to_reference_structure(mapping_df, reference_struct):
    '''IN: pandas dataframe of all structures mapped to alingment position,
           reference structure (included in the dataframe)
    OUT: same dataframe with added columns: reference structure position'''
    ref_df = mapping_df[mapping_df.pdb_id == reference_struct][['pdb', 'alignment_pos']]
    ref_df.columns = ['ref_pdb', 'alignment_pos']
    mapping_df = pd.merge(mapping_df, ref_df, on='alignment_pos', how='left')
    mapping_df = mapping_df.sort_index(by=['pdb_id', 'resnum'])
    mapping_df = mapping_df.reset_index(drop=True)
    return mapping_df


mapping = pd.DataFrame(columns = ['pdb_id', 'resnum', 'pdb', 'alignment_pos',
                                  'aa'])
for filename in os.listdir(args.processed_pdb_dir):
    print filename
    pdbID = filename.split('.')[0]
    parser = PDB.PDBParser()
    struct = parser.get_structure(pdbID, os.path.join(args.processed_pdb_dir,
                                                      filename))
    res_list = list(struct.get_residues())
    alignment_ref_list = map_to_alignment(args.reference_alignment, filename)
    for i in range(0, len(res_list)):
        mapping = mapping.append({'pdb_id': pdbID,
                                  'resnum': (i + 1),
                                  'pdb': res_list[i].id[1],
                                  'alignment_pos': alignment_ref_list[i],
                                  'aa': res_list[i].resname},
                                  ignore_index=True)

mapping = map_to_reference_structure(mapping, args.reference_structure)

# WRITE MAPPING FILE
mapping.to_csv('results/mapping.csv', index=False)

