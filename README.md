---
output: html_document
---

# Consensus residue contact calculator


## What does this software do?

The program calculates residue contact networks for a number of related protein structures from the protein data base (PDB) provided by the user. Furthermore, for every residue contact, the software determines which fraction of structures have an equivalent contact. In order to identify common residue contacts, the protein structures provided by the user have to be related (The software will finish either way but won't find any common residue contacts if unrelated protein structures are provided). In order to identify structurally equivalent residues, the user has to provide an alignment of the sequences of the protein structures.

The user can chose a particular structure of interest. The results show, for every residue contact in the structure, the fraction of related protein structures which have an equivalent contact. This information can be used to gauge the structural importance of residue contacts thus aiding interpretation of protein structures. For instance, a residue contact present in all members of a given protein family can be expected to be important for the common fold of this protein family. Moreover, the results allow to distinguish between features unique to a given protein structure of interest and features shared across the whole protein family.

Finally, the software automatically generates a report which visualizes the results and can serve as a base for further analysis.

For further information about the method of consensus residue contact networks, please refer to the documentation and the papers cited therein.  **insert reference(s)**





## Software requirements

The program was developed and has been tested on Linux. It should work on any Linux-distribution.

Before running the program, please make sure you have the following dependencies installed on your computer:

* Python version 2.7 or newer
* IPython
* Python modules: os, sys, argparse, numpy, pandas, networkx, Bio (specifically SeqIO and PDB)
* R version 3.2.5 or newer
* R packages: plyr, ggplot2, Biostrings


## Input data

For simplicity, it is recommended to create a dircetory 'data' in the directory of the software (Consensus-residue-contact-calculator-x.x). As an example, see the directory 'test_data'. Make sure your data-directory has the following content:

* PDB-files: Select protein structures of interest from a given protein family and download the corresponding PDB-files. Put them in a directory which contains no other files.
* reference alignment: Create an alignment of all sequences of the structures to be analyised. Make sure the names in the alignment are the same as the names in the names of the PDB-files described above. This alignment is used by the software to identify structurally equivalent residues. Therefore, for proteins with low sequence similarity, it is recommended to use structural alignment (and get the sequence alignment form the structural alignment) rather than sequence alignment. However, most protein families have sufficient sequence similarity to use 
* pdb_chain_pfam.csv: You can simply use the file from the directory 'test_data' or download the latest version from SIFTS-database: from https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html


## How do I get set up?

### Summary

1) Download the repository as zip-file from: https://github.com/raph333/Consensus-residue-contact-calculator/releases/
2) Unpack the zip-file:  
unzip Consensus-residue-contact-calculator-x.x  # x.x being the current release number
3) Test with example data:  
cd Consensus-residue-contact-calculator-x.x  # x.x being the current release number  
bash runall.sh test_data/raw_pdb_files test_data/ras_reference_alignment.fa 1g16 PF00071
4) If the script has finished, inspect the folder 'results' to see what to expect from the analysis. In case of an error, please address follow the instructions provided by the error-message (e.g. install required software).
5) Start the script 'runall' with your own data. Run  
bash runall.sh --help  
or read the next section for more information about the required arguments.


### Usage



## Who do I talk to?

For questions regarding the usage of the software, bug-reports or suggestions, please don't hesitate write an e-mail to raphae1peer@gmail.com.
