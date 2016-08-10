'''
------------------------------------------------------------------------------
AUTHOR: Raphael Peer, raphael1peer@gmail.com

PURPOSE:
1) The script tries to import all python modules required by the python scripts
of the 'Consensus residuce contact calculator' and exits if one of the modules
is not installed.
2) The script checks the validity of the input data provided by the user
a) PDB-files
b) reference alignment (in fasta format)
c) csv-file from the SIFTS-database (which contains a list of pdb-IDs, chains
and the respective Pfam domain)
If a problem with the input data is encountered, the an error message is given
and the script and exits.
------------------------------------------------------------------------------
'''

# Check if all python modules required for the pipeline are installed,
# otherwise exit script:
print('Checking software requirements:')
try:
    import os, sys, argparse, numpy, pandas, networkx
    from Bio import SeqIO, PDB
    print('All required python modules are installed.')
except ImportError:
    print('One of the python modules required for this software is missing.')
    print('Please installe the missing module:')
    import os, sys, argparse, numpy, pandas, networkx
    from Bio import SeqIO, PDB


parser = argparse.ArgumentParser()
parser.add_argument('raw_pdb_dir', help='Directory with PDB-structures '
                    'for calculation of residue contact network')
parser.add_argument('reference_alignment', help='Alignment of the sequences '
                    'of the structures for which residue contact networks '
                    'were created. See the documentation for information '
                    'about the requirements of such an alignment')
parser.add_argument('sifts_chain_pfam', help='SIFTS-file "pdb_chain_pfam.csv" for '
                    'finding chains withing the PDB-structures which contain '
                    'the Pfam-domain of interest')
try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(1)


def check_pdb_files(pdb_file_dir):
    '''Try to parse all PDB-files in a directory with Bio.PDB;
    Exit script if parsing of fails.'''
    print('\nChecking PDB-files in directory "%s":' % pdb_file_dir)
    filecount = 0
    for filename in os.listdir(pdb_file_dir):
        if not filename.endswith('.pdb'):
            print('Unexpected filenname: %s.' % filename)
            print('Please make sure, that all PDB-files end with ".pdb".')
            sys.exit(1)
        pdb_id = filename.split('.')[0]
        struct = PDB.PDBParser().get_structure(pdb_id,
                                        os.path.join(pdb_file_dir, filename))
        #print pdb_id, len(list(struct.get_chains()))
        if len(list(struct.get_chains())) < 1:
            print('Problem encountered in file "%s": not a valid '
                'PDB-file.\nPlease remove it from the directory.' % filename)
            sys.exit(1)
        filecount += 1
    print('%s valid PDB-files found.' % filecount)


def check_SIFTS_file(sifts_pdb_pfam_file):
    '''Try to read csv-file with pandas; Exist script if reading fails'''
    print('\nChecking SIFTS-file: %s' % sifts_pdb_pfam_file)
    print('(required to identify PDB-chains with Pfam-domain of interest)')
    try:
        chain_pfam = pandas.read_csv(sifts_pdb_pfam_file, comment='#')
        print('File present: %s' % sifts_pdb_pfam_file)
        return chain_pfam
    except:
        print('File %s not present or corrupted.' % sifts_pdb_pfam_file)
        print('Please check the file. If necessary, you can download the '
            'correct file here: '
            'https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html')
        sys.exit(1)


def check_reference_alignment(reference_alignment, pdb_files_dir):
    '''Try to read alignment (with Bio.SeqIO); Exit if reading fails;
    Print warining message (but do not exit) if a structure is not in the
    alignment'''
    print('\nChecking reference alignment: %s' % reference_alignment)
    try:
        ref_alignment = SeqIO.parse(reference_alignment, 'fasta')
        print('Valid alignment present: %s' % reference_alignment)
    except:
        print('Alignment not present or not in required format (fasta).')
        print('Please check the file "%s".' % reference_alignment)
        sys.exit(1)
    
    seq_names_alignment = SeqIO.to_dict(ref_alignment).keys()
    structures = os.listdir(pdb_files_dir)
    not_in_alignment = [x for x in structures if x not in seq_names_alignment]
    if len(not_in_alignment) > 0:
        print('WARNING: The sequence of the following %s PDB-structures is no '
            'present in the reference alignment and will be excluded from '
            'analysis:\n%s' % (len(not_in_alignment), not_in_alignment))
        print('If the dataset is large enought without these structures, this '
            'is not a problem.')
        print('Note: The sequences in the alignment must have the same names '
            'as their corresponding structures in the directory "%s".'
            % pdb_files_dir)
    else:
        print('All structures are present in the reference alignment.')


# CHECK INPUT DATA:
check_pdb_files(args.raw_pdb_dir)
check_SIFTS_file(args.sifts_chain_pfam)
check_reference_alignment(args.reference_alignment, args.raw_pdb_dir)
