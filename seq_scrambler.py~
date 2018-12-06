""" Assemble a MSA to a especific position within coordinates given by a matrix. """

# Libraries:

import os
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from collections import OrderedDict as odict

# Constants:

WD = os.chdir('/home/earaujo/Repositories/HLA/workbench/gen_a/')
MSA = SeqIO.parse('CD8.last.fasta', 'fasta')
GENOME = np.load('genomes/genome.7000.npy')

# Functions:

record_dict = odict(); count = 0
for record in MSA:
    record_dict[str(count)] = record
    count += 1
    
to_write = []
for index in range(0, len(GENOME)):
    to_write.append(record_dict[str(GENOME[index])])
    
SeqIO.write(to_write, 'result.fasta', 'fasta')



