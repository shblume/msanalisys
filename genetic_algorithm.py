import os
import numpy as np
import random as rd
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dtime
from Bio import AlignIO
from scipy import spatial
from joblib import Parallel, delayed
from functions import *

# USER DEFINED #############################################
GENERATIONS = 10000
RESTART = 0
SAVEFREQ = 1
LAMBDA = 0.001
OFFSET = 95
theta = 1.0 # Between 0.0 and 1.0
q = 21

# Parametros de paralelização
MULT_CHAIN = 0
NUM_CF = 25
NUM_CM = 10
NUM_CQ = 5
MUT_INIT = 2
MUT_CF = 2
MUT_CM = 4
MUT_CQ = 50
RAND = 1
NUM_CORES = 4
CRITERIA = 0
MINIMIZATION = 0
OPTMIZE = 2 # 0: MI - 1: MI_PP_by_Meff - 2: H(x,y)/SP

# Files
GENOMES_PATH = "genomes"
BEST_PATH =     "best.mi"
RESTART_FILE = "{}/genome.{}.npy".format(GENOMES_PATH, RESTART)
SYSTEM = 'HLA_CD8'
MSA_TOX = "HLA_GA_edited.fa"
MSA_NAV = "CD8.last.fa"

###########################################################

# CREATING DIRECTORIES
if not os.path.exists(GENOMES_PATH):
    os.makedirs(GENOMES_PATH)

#######################################################################################################
#######################################################################################################

# LOADING FILES
ics_a = [1,2,4,5,6,9,10,13,16,19,20,33,48,57,61,62,63,67,70,75,77,79]
ics_b = list(range(OFFSET,OFFSET + 108))
ics = np.concatenate((ics_a, ics_b))

#READ FASTA-FORMAT MSA INTO BIOPYTHON ALIGNMENT OBJECT
tox_handle = open(MSA_TOX, "r")
nav_handle = open(MSA_NAV, "r")
msa_tox = AlignIO.read(tox_handle, "fasta")
msa_nav = AlignIO.read(nav_handle, "fasta")
seqnumber = len(msa_tox)
seqlength = len(msa_tox[0]) + len(msa_nav[0])

# ENCODE MSA
aminos = {"A": 0, "R": 1, "N": 2, "D": 3, "Q": 4,
          "E": 5, "G": 6, "H": 7, "L": 8, "K": 9,
          "M":10, "F":11, "S":12, "T":13, "W":14,
          "Y":15, "C":16, "I":17, "P":18, "V":19,
          "-":20, ".":20, "B": 2, "Z": 4, "X":20, "J":20}

encoded_msa0=np.empty((seqnumber,seqlength),dtype=int)
for (x,i), A in np.ndenumerate(msa_tox):
    encoded_msa0[x,i]=aminos[A.upper()]
for (x,i), A in np.ndenumerate(msa_nav):
    encoded_msa0[x,i+OFFSET]=aminos[A.upper()]

seqs_b = encoded_msa0[:,OFFSET:]

# defining sizes
nA = len(ics_a)
nB = len(ics_b)
seqlength = nA+nB

# defining contacts
pairs = []
idx_pairs = []
for i, ai in enumerate(ics_a):
    for j, aj in enumerate(ics_b):
        pairs.append((ai,aj))
        idx_pairs.append((i,j))

nP = len(pairs) # important for SP function

# BUILD GENOMES LIST
genomes = []

# READING GENOME TO START SIMULATION
if RESTART == 0:
    for j in range(NUM_CORES):
        g = list(range(seqnumber))
        rd.shuffle(g)
        genomes.append(g)
    START = 0
else:
    START = RESTART + 1
    handle = np.load(RESTART_FILE)
    for j in range(NUM_CORES):
        genomes.append(list(handle))

############################################################
###########################################################
 
def Fitness(g):
    encoded_msa, Meff = Codemsa(g, encoded_msa0, seqs_b, OFFSET, theta)
    sitefreq = Sitefreq(encoded_msa, Meff, ics, nA+nB, q, LAMBDA)
    pairfreq = Pairfreq(encoded_msa, Meff, ics, nA+nB, sitefreq, q, LAMBDA)
    mi, h = information(sitefreq, pairfreq, nP, pairs, idx_pairs, q)

    fitness = np.sum(mi)

    return [g, fitness]
   

################################################################
###############################################################

# GA loop 
#dataframe = open('HLA.genome.result.dat', 'a')
genomes = genomes[1:].copy()
for generation in range(START, GENERATIONS):
    # Restarting the mutation rate
    MUTATIONS = MUT_INIT
    #   MUTATIONS  #
    for n in range(NUM_CORES - 1):
        indexes = []
        
        # Defining the best genome to mutate
        if MINIMIZATION == 1:
            gen = genomes[0].copy()
        else:
            gen = genomes[-1].copy()

        #   MUTATION  #
        for j in range(MUTATIONS):
            a, b = rd.sample(gen,2)
            while (a,b) in indexes:
                a, b = rd.sample(gen,2)

            gen[int(a)], gen[int(b)] = gen[int(b)], gen[int(a)] # changing indices of genome
            indexes.append((a,b))
            indexes.append((b,a))

        # Transfering gen with mutation to population (genomes)
        if MINIMIZATION == 1:
            genomes[n+1] = gen.copy()
        else:
            genomes[n] = gen.copy()

        # changin mutation rate
        if MULT_CHAIN == 1:
            if n < NUM_CF:
                MUTATIONS = MUT_CF
            if n >= NUM_CF and n < NUM_CF + NUM_CM:
                MUTATIONS = MUT_CM
            if n >= NUM_CF + NUM_CM:
                MUTATIONS = MUT_CQ
        else:
            MUTATIONS += 0
    
    # CALCULATING THE OPTIMIZATION FUNCTION WITH MULTI THREADS
    results = Parallel(n_jobs=NUM_CORES)(delayed(Fitness)(g) for g in genomes)

    genomes_pre = genomes.copy()
   
    # DELTA
    results_selected = []
    if generation != 0 and CRITERIA == 1:
        for k in results:
            if best_pre[2] != k[2] and best_pre[3] != k[3]:
                delta_frac = (k[2] - best_pre[2]) / (best_pre[3] - k[3])
                if delta_frac > HSP_FRAC:
                    results_selected.append(k)
            else:
                results_selected.append(k)
    else:
        results_selected = results

    genomes = []
    results_sorted = sorted(results_selected, key=lambda results : results[1])
    for result in results_sorted:
        genomes.append(result[0])
    best_pre = results_sorted[0]

    while len(genomes) < NUM_CORES:
        genomes.append(genomes[0])

    # SAVING DATA
    if MINIMIZATION == 1:
        
        if generation % SAVEFREQ == 0:
            np.save("{}/genome.{}".format(GENOMES_PATH, generation), genomes[0])
    else:
        if generation % SAVEFREQ == 0:
            np.save("{}/genome.{}".format(GENOMES_PATH, generation), genomes[-1])

    print("Generation {}".format(generation))
#    dataframe.write('\n GEN {}:         '.format(generation))
    for k in results_sorted:
        print(k[1])
        np.save('{}/best.{}'.format(BEST_PATH, generation), np.array(k[1]))
#        dataframe.write('{}; {}; {}'.format(k[1][0],k[1][1],k[1][2]))
                
