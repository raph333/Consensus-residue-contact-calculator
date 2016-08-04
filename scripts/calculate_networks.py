import argparse
import sys
import os
import numpy as np
import pandas as pd
import networkx as nx
from Bio import PDB

#import time
#start =  time.time()

parser = argparse.ArgumentParser()
parser.add_argument('processed_pdb_dir', help='directory with processed '
                    'PDB-structures for calculation of residue contact network')
parser.add_argument('cutoff', nargs='?',type=int, const=5, default=5,
                    help='Any two residues which have at least one pair of '
                    'atoms within this distance are considered to make a '
                    'contact. If no argument is provided, the default value '
                    'of 5 Angstrom is used.')
try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(1)
print 'atomic distance cutoff: %s Angstrom' % args.cutoff


def min_atomic_distance(resA, resB):
    '''IN: two residues (Bio.PDB objects)
    OUT: shortest distance between any two atoms'''
    atoms1 = list(resA.get_iterator())
    atoms2 = list(resB.get_iterator())
    #print "\n\nresidue 1: %s %s: %s" % (resA.id[1], resA.resname, atoms1)
    #print "residue 2: %s %s: %s" % (resB.id[1], resB.resname, atoms2)
    atomic_dist_matrix = np.zeros((len(atoms1), len(atoms2)), dtype=('f4'))
    for i in range(0, len(atoms1)):
        for j in range(0, len(atoms2)):
            atomic_dist_matrix[i,j] = atoms1[i] - atoms2[j]
    #print "%s" % atomic_dist_matrix
    return atomic_dist_matrix.min()  # minimal atomic distance


def residue_distance_matrix(directory, filename):
    '''IN: input directory and pdb-filename
    OUT: 2D numpy array of inter-residue distances (NA = large distance)'''        
    struct = PDB.PDBParser().get_structure(filename, os.path.join(directory, filename))
    res_list = list(struct.get_residues())
    dist_matrix = np.zeros((len(res_list), len(res_list)), dtype=('f4'))
    dist_matrix[:] = np.nan
    for i in range(0, len(res_list)):
        for j in range(i, len(res_list)):  # start at i: only half matrix
            # if the CA-CA distance is above 15A, don't calculate
            # any to any atom distances: reduces runtime to 1/3
            try:
                CA_dist = res_list[i]['CA'] - res_list[j]['CA']
                if CA_dist <= 15:
                    atomic_distance = min_atomic_distance(res_list[i], res_list[j])
                    dist_matrix[i,j] = atomic_distance
            except:  # if residue does not have CA-atom coordinates
                pass
    return(dist_matrix)


def matrix_to_dataframe(matrix, pdb_id):
    '''IN: all against all residue matrix ('1' means contact)
    OUT: dataframe with one contact per line (residue A, residue B)'''
    nw = nx.from_numpy_matrix(contact_matrix)
    # network has sequential residue numbering starting at 0 (no pdb-numbers)
    df = pd.DataFrame(nw.edges())  # only edges with weight '1' are taken
    df['pdb_id'] = pdb_id
    df.columns = ['res_A', 'res_B', 'pdb_id']
    df.res_A += 1  # residue numbering should start at 1
    df.res_B += 1
    return(df)


filecounter = 0
for filename in os.listdir(args.processed_pdb_dir):
#    if not filename.endswith('.pdb'):
#        continue  # ignore non-PDB files
    pdb_id = filename.split('.')[0]
    distance_matrix = residue_distance_matrix(args.processed_pdb_dir, filename)
    contact_matrix = (distance_matrix < args.cutoff) & (distance_matrix > 0).astype(int)
    nw = matrix_to_dataframe(contact_matrix, pdb_id)
    filecounter += 1
    if filecounter == 1:
        networks = nw  # initizalize big dataframe to store all networks
    else:
        networks = pd.concat([networks, nw])
    print '%s: %s' % (filecounter, pdb_id)


# WRITE ALL NETWORKS TO FILE
networks = networks[['pdb_id', 'res_A', 'res_B']]
networks = networks.sort_index(by=['pdb_id', 'res_A', 'res_B'])
networks.to_csv('results/raw_networks.csv', index=False)

#end =  time.time()
print "Number of residue contact networks calcuated: %s" % filecounter
#print "Runtime: %s minutes" % round(((end - start) / float(60)), 1)

