""" This program identifies and enumerates sequences from a clustal file about another fasta archive. """

# Libraries imported:
import Bio
from Bio import SeqIO
from Bio import AlignIO

# Functions defined:

def enumerator(align,ref_list,number):
    to_be_written = []; number = str(number); counter = 0
    for arec in AlignIO.read(align+'.aln','clustal'):
        aln_id = arec.id
        for srec in SeqIO.parse(ref_list+'.fa','fasta'):
            seq_id = srec.id
            counter += 1
            if counter // 1000 == 0:
                print(' Please, wait ...')
            if arec.seq[-1] == number or arec.seq[-1] == number:
                continue
            else:
                if aln_id.startswith(seq_id):
                    to_be_written.append(arec + number)
    
    return to_be_written

# Execution:

aln_name = input(' Enter the alignment to be enumerated: ')
references = []; counter = 1 
ref = input(' Enter the reference for the number 0: '); references.append(ref)
while ref != 'done' or ref != 'Done' or ref != 'DONE':
    ref = input(' Enter the reference for the number {}: '.format(counter))
    if ref == 'done' or ref == 'Done' or ref == 'DONE':
        break
    else:
        counter += 1
        references.append(ref)
print(references) 
gen_list = []
for index in range(0, len(references)):
    gen_list += enumerator(aln_name, references[index], index)

SeqIO.write(gen_list,aln_name+'.enum.aln','clustal')

print(' Program succesful! Check your file!.')

print(' Program closing.')
