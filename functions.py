import numpy as np
from Bio import AlignIO
from scipy import spatial
USE_W = 0
if USE_W:
    W = np.loadtxt("wfunc2.dat")
else:
    W = np.ones((21,21))

def Corr(g, distances_A, distances_B):
    x = distances_A
    y = distances_B
    num = 0
    den_r = 0
    den_s = 0
    R_av = np.average(distances_A)
    S_av = np.average(distances_B)
    for i, i_b in enumerate(g):
        for j, j_b in enumerate(g):
            num += (x[i,j] - R_av) * (y[int(i_b),int(j_b)] - S_av)
            den_r += (x[i,j] - R_av)**2
            den_s += (y[int(i_b),int(j_b)] - S_av)**2
    r = num / (np.sqrt(den_r) * np.sqrt(den_s))
    return(r)

def Codemsa(g, encoded_msa0, seqs_b, OFFSET, theta):
    g = np.array(g).astype(int)
    encoded_msa = encoded_msa0.copy()
    encoded_msa[:,OFFSET:] = seqs_b[g]

    #WEIGHT SEQUENCES
    hammdist = spatial.distance.pdist(encoded_msa, 'hamming')
    weight_matrix = spatial.distance.squareform(hammdist < (1.0- theta))
    weight = 1.0 / (np.sum(weight_matrix, axis = 1) + 1.0)
    Meff = np.sum(weight)
    return encoded_msa, Meff

def Sitefreq(encoded_msa, Meff, ics, nA, q, LAMBDA):
    sitefreq = np.empty((nA,q),dtype=float)
    for i, col_i in enumerate(ics):
        for aa in range(q):
            sitefreq[i,aa] = np.sum(np.equal(encoded_msa[:,col_i],aa))/Meff
    sitefreq = (1-LAMBDA)*sitefreq + LAMBDA/q
    return sitefreq

def Entropy(sitefreq,nA,nB):
    ent = np.zeros((nA+nB),dtype=float)
    for i in range(nA+nB):
        ent[i] = -np.sum(sitefreq[i,:]*np.log(sitefreq[i,:]))
    return ent

def cantor(x, y):
    return (x + y) * (x + y + 1) / 2 + y

def Pairfreq(encoded_msa, Meff, ics, nP, sitefreq, q, LAMBDA):
    pairfreq = np.zeros((nP,q,nP,q),dtype=float)
    for i, col_i in enumerate(ics):
        for j, col_j in enumerate(ics):
            c = cantor(encoded_msa[:,col_i],encoded_msa[:,col_j])
            unique,aaIdx = np.unique(c,True)
            for x,item in enumerate(unique):
                pairfreq[i,encoded_msa[aaIdx[x],col_i],j,encoded_msa[aaIdx[x],col_j]] = np.sum(np.equal(c,item))

    pairfreq /= Meff
    pairfreq = (1-LAMBDA)*pairfreq + LAMBDA/(q*q)

    for i in range(nP):
        for am_i in range(q):
            for am_j in range(q):
                if (am_i==am_j):
                    pairfreq[i,am_i,i,am_j] = sitefreq[i,am_i]
                else:
                    pairfreq[i,am_i,i,am_j] = 0.0
    return pairfreq

def information(sitefreq, pairfreq, nP, pairs, idx_pairs, q):
    H_xy_matrix_pp = np.empty((nP),dtype=float)
    mi_matrix_pp = np.empty((nP),dtype=float)
    tiny = 1e-10
    for k,(col_i,col_j) in enumerate(pairs):
        i = idx_pairs[k][0]
        j = idx_pairs[k][1]
        pij = pairfreq[i,:,j,:]
        pi = np.transpose(np.broadcast_to(sitefreq[i,:],(q,q)))
        pj = np.broadcast_to(sitefreq[j,:],(q,q))
        mi_matrix_pp[k] = np.sum(pij*np.log(pij/(pi*pj) + tiny))
        H_xy_matrix_pp[k] = -np.sum(pij*np.log(pij + tiny))

    return mi_matrix_pp, np.sum(H_xy_matrix_pp)

def calc_delta_f(pairfreq1, pairfreq2, ics, pairs):
    for i,col_i in enumerate(ics):
        for j,col_j in enumerate(ics):
            if (col_i,col_j) not in pairs:
                pairfreq1[i,:,j,:] = 0
                pairfreq2[i,:,j,:] = 0
    delta_f = pairfreq1 - pairfreq2
    delta_f *= delta_f
    delta_f = np.sqrt(np.sum(delta_f))

    return delta_f
