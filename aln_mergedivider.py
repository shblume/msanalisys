""" This program divides the MSA based on the last column, generated by the 'aln_seq_identifier.py' program or merges two clustal files independently """

# Libraries imported:

from Bio import AlignIO
from Bio import SeqIO

# Functions defined: 

def divider(align,a_name,b_name):
    a_list = []; b_list =[]
    for record in AlignIO.read(align+'.aln','clustal'):
        if record.seq[-1] == '0':
            a_list.append(record)
        elif record.seq[-1] == '1':
            b_list.append(record)
    SeqIO.write(a_list,a_name+'.aln','clustal')
    SeqIO.write(b_list,b_name+'.aln','clustal')

def merger(a_align,b_align,name):
    wanted = []
    for a_record in AlignIO.read(a_align+'.aln','clustal'):
        wanted.append(a_record)
    for b_record in AlignIO.read(b_align+'.aln','clustal'):
        wanted.append(b_record)
    SeqIO.write(wanted,name+'.aln','clustal')

# Execution:

startpoint = input(' Do you want to merge (m) or divide (d)? [m/d]: ')
if startpoint == 'd' or startpoint == 'D':
    aln = input(' Enter the clustal file: ')
    first_name = input(' Enter the first file name: ')
    second_name = input(' Enter the second file name: ')
    divider(aln,first_name,second_name)
    print(' Program successful! Check your file.')
elif startpoint == 'm' or startpoint == 'M':
    first_aln = input(' Enter the first clustal file: ')
    second_aln = input(' Enter the second clustal file: ')
    name = input(' Enter the new alignment name: ')
    merger(first_aln,second_aln,name)
    print(' Program successsful! Check your file.')
print(' Program closing.')
