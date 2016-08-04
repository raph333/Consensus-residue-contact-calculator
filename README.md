---
output: pdf_document
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


## How do I get set up?

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
