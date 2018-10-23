""" This program removes from a fasta archive the unwanted sequences. """

# Libraries imported:

import re
from Bio import SeqIO

# Functions defined:

def name_remover(fname, ref):
    counter = 0
    wanted_seq = []
    for record in SeqIO.parse(fname+'.fa', 'fasta'):
        identif = record.id; #print(identif)
        if identif.startswith(ref): 
            counter += 1
        else:
            wanted_seq.append(record); 
    SeqIO.write(wanted_seq, fname+'.fa','fasta')            
    return counter
    
def size_remover(fname, size=10, par='less'):
    wanted_seq = []; counter = 0
    for record in SeqIO.parse(fname+'.fa', 'fasta'):
        if par == 'less':
            if len(record.seq) <= size:
                counter += 1
                continue
            else:
                wanted_seq.append(record)
        if par == 'more':
            if len(record.seq) >= size:
                counter += 1
                continue
            else:
                wanted_seq.append(record)
    SeqIO.write(wanted_seq, fname+'.fa', 'fasta')
    return counter

# Execution:

start = input(' Do you want to use a archive as reference? [y/n]: ')
archive = input(' Please, insert the base fasta file name: ')

if start == 'y' or start == 'Y':
    reference = input (' Please, insert the reference fasta file: ')
    for rem_record in SeqIO.parse(reference+'.fa','fasta'):
        removed = remover(archive,rem_record.id)
        if removed == 0:
            print(' No sequence with the reference archive "{}.fa" were found in the "{}.fa" archive.'.format(reference,rem_record.id))
        else:
            print(' There were {} removed sequences with the "{}" reference.'.format(removed,rem_record.id))
    print('\n Program successful! Check your file.')
       
else:
    reference = None
    while reference != 'done' or reference != 'Done' or reference != 'DONE':
        reference = input(' Please, insert the exclusion reference: ')
        if reference == 'done' or reference == 'Done' or reference == 'DONE':
            print(' Program successful! Check your file.')
            break
        elif reference[0] == '-':
            if reference[1:5] == 'less':
                number = reference[6:-1]+reference[-1]; number = int(number)
                removed = size_remover(archive, size=number)
                if removed == 0:
                    print('\n No sequence with less than {} characters were found.'.format(number))
                else:
                    print('\n {} sequences with less than {} characters were removed.'.format(removed,number))
            elif reference[1:5] == 'more':
                number = reference[6:-1]+reference[-1]; number = int(number)
                removed = size_remover(archive, size=number, par='more')
                if removed == 0:
                    print('\n No sequence with more than {} characters were found.'.format(number))
                else:
                    print('\n {} sequences with more than {} characters were found.'.format(removed,number))
            elif reference[1:5] == 'diff':
                archive = reference[6:-1]+reference-1
        else:
            removed = name_remover(archive,reference)
            if removed == 0:
                print('\n No sequence with the reference "{}" were found in the {}.fa arquive.\n'.format(reference,archive))
            else:
                print('\n There were {} removed sequences with the "{}" reference.\n'.format(removed,reference))

print(' Program closing.')
