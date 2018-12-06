""" This program identifies all identical sequences between two fasta files. """

# Libraries imported:

from Bio import SeqIO

# Functions defined:

def spotter(a_file,b_file):
    repeated = [];
    for a_record in SeqIO.parse(a_file, 'fasta'):
        a_ident = a_record.id
        for b_record in SeqIO.parse(b_file, 'fasta'):
            b_ident = b_record.id
            if a_ident == b_ident:
                repeated.append(a_record)
    SeqIO.write(repeated, 'repeated_sequences.fa','fasta')
    return repeated

# Execution:

a_inp = input(' Please, enter the first .fasta archive: '); a_fasta = a_inp+'.fa'
b_inp = input(' Enter another .fasta archive: '); b_fasta = b_inp+'.fa'
rep_list = spotter(a_fasta,b_fasta)
if rep_list == []:
    print('\n There were no repeated sequences on the two given fasta files.')

else:
    print('\n There were {} repeated sequences:'.format(len(rep_list)))
