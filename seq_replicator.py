''' Replicates a given sequence or a random number of sequences. '''

# Libraries imported:

import numpy as np
import random as rd
import Bio
from Bio import SeqIO
from Bio import AlignIO

# Function defined:

def replicator(aln, n=17):
    id_list = []; rec_list = []; choosen = [] ; index = 0
    for record in SeqIO.parse(aln+'.fa', 'fasta'):
        id_list.append(record.id)
    while index < n:
        choosen.append(rd.choice(id_list))
        index += 1
    for record in SeqIO.parse(aln+'.fa', 'fasta'):
        for item in choosen:
            if record.id == item:
                rec_list.append(record); rec_list.append(record)
            else:
                rec_list.append(record)

    return rec_list
    
# Program execution:

fil_name = input(' Insert the fasta: ')
number = input(' Insert how many replications you want: '); number = int(number)

w_list = replicator(fil_name, number)
SeqIO.write(w_list, fil_name+'.dpc.fa', 'fasta')

print('\n Program sucessful!')
