
# coding: utf-8

# KDD’16 "Joint Community and Structural Hole Spanner Detection via Harmonic Modularity” @Lifang 07/16/2016
import sys
import numpy as np
import json, os
import scipy.optimize as spo
import scipy.sparse as sps
import scipy.linalg as spl
import networkx as nx
import random
import scipy.io as sio
from sklearn.cluster import KMeans
from sklearn import metrics
import scipy.stats as stat
from sklearn.preprocessing import normalize
from scipy.cluster.vq import kmeans, vq, kmeans2
from copy import deepcopy
from collections import Counter
import matplotlib as mpl
import matplotlib.pyplot as plt
# %matplotlib inline
eps=2.220446049250313e-16

# Generate an orthogonal matrix with a given matrix
def sym(w):
    return w.dot(spl.inv(spl.sqrtm(w.T.dot(w))))

def load_adj_matrix(fname):
    adj_tuples =  np.loadtxt(fname,dtype='int')
    n = adj_tuples[0,0] # number of nodes
    vals = np.array([1] * len(adj_tuples[1:])) # values = 1
    max_id = max(max(adj_tuples[1:,0]),max(adj_tuples[1:,1])) + 1
    A = sps.csr_matrix((vals,(adj_tuples[1:,0], adj_tuples[1:,1])), shape=(max_id,max_id))
    A = A + A.T
    return sps.csr_matrix(A)

def load_labels(fname):
    fin=open(fname,'r')
    label = 0
    gdclus = {}
    cluster_dict = {}
    for line in fin:
        nodes = line.split('\t')
        for node in nodes:
            if int(node) in gdclus:
                gdclus[int(node)].append(label)
            else:
                gdclus[int(node)] = [label]
        cluster_dict[label] = [int(node) for node in nodes]
        label += 1
    return gdclus, cluster_dict

# classifify SHS using majority voting
def majority_voting(votes):
    C = Counter(votes)
    pairs = C.most_common(2)
    if len(pairs)==0:
        return 0
    if pairs[0][0] > 0:
        return pairs[0][0]
    elif len(pairs)>1:
        return pairs[1][0]
    else:
        return 0

def label_by_neighbors(AdjMat,labels):
    assert (AdjMat.shape[0] == len(labels)), "dimensions are not equal"
#     print labels
    unlabeled_idx = (labels==0)
    num_unlabeled = sum(unlabeled_idx)
    count = 0
    while num_unlabeled > 0:
        idxs = np.array(np.nonzero(unlabeled_idx)[0])
        next_labels = np.zeros(len(labels))        
        for idx in idxs:
            neighbors = np.nonzero(AdjMat[idx,:] > 0)[1]
            if len(neighbors)==0:
                next_labels[idx] = majority_voting(labels)
            neighbor_labels = labels[neighbors]
#             print idx, neighbors, neighbor_labels
            next_labels[idx] = majority_voting(neighbor_labels)
        labels[idxs] = next_labels[idxs]
        unlabeled_idx = (labels==0)
        num_unlabeled = sum(unlabeled_idx)
#         print num_unlabeled
    return labels

def save_labels(fname, labels):
    fout = open(fname,'w')
    unique_labels = np.unique(labels)
    for community in unique_labels:
        nodes = np.nonzero(labels == community)[0]
        #     change label starting index from 0 to 1 (as in mtx file format)
        community_str = '\t'.join([str(node+1) for node in nodes])
        if community == unique_labels[-1]:
            fout.write(community_str)
        else:
            fout.write(community_str+'\n')
    fout.close()

import scipy.stats as stat
def avg_entropy(predicted_labels, actual_labels):
    actual_labels_dict = {}
    predicted_labels_dict = {}
    for label in np.unique(actual_labels):
        actual_labels_dict[label] = np.nonzero(actual_labels==label)[0]
    for label in np.unique(predicted_labels):
        predicted_labels_dict[label] = np.nonzero(predicted_labels==label)[0]
    avg_value = 0
    N = len(predicted_labels)
    # store entropy for each community
    for label, items in predicted_labels_dict.iteritems():
        N_i = float(len(items))
        p_i = []
        for label2, items2  in actual_labels_dict.iteritems():
            common = set(items.tolist()).intersection(set(items2.tolist()))
            p_ij = float(len(common))/ N_i
            p_i.append(p_ij)
        entropy_i = stat.entropy(p_i)
        avg_value += entropy_i * (N_i / float(N))
    return avg_value

## Loading data ......

datapath= os.getcwd() + '/data/'
adj_name = 'karatern'     # Note: the data index starts from 0, not from 1
community_name = 'karatecrn'

#fname = datapath+'karatern.txt'
fname = datapath + adj_name +'.txt'
A_mat = load_adj_matrix(fname)

#fname = datapath+'karatecrn.txt'
fname = datapath + community_name +'.txt'
gdclus, cluster_dict = load_labels(fname)
ground_truth_labels = list(gdclus.values())

## Algorithm: Harmonic Modularity (HAM)

# Initialize
A = A_mat        # adjacency matrix
n = A.shape[0]   # the number of nodes
c=len(np.unique(ground_truth_labels))  # the number of communities
topk=3           # the number of selected SH spanner
epsilon=1e-4       # smoothing value: epsilon
max_iter = 50    # maximum iteration value
seeeed = 5433
np.random.seed(seeeed)

invD=sps.diags(np.array(A.sum(axis=0))[0,:]**(-1.0),0)   # degree matrix D^-1
L = (invD.dot(A)-sps.identity(n)).tocsr()     # Laplacian matrix L= D^-1 * A - I

F=sym(np.random.random((n,c)))               # Generate a random orthogonal matrix

# Main Body
for step in range(max_iter):
    Q=sps.identity(n).tocsr()
    P = L.dot(F)
    for i in range(n):
        Q[i,i]=0.5/(spl.norm(P[i,:])+epsilon)
    
    R=L.T.dot(Q).dot(L)
    
    W,V=np.linalg.eigh(R.todense())
    Wsort=np.argsort(W)      # sort from smallest to largest   
    F=V[:,Wsort[0:c]]        # select the smallest eigenvectors

# find SH spanner    
SH=np.zeros((n,))
for i in range(n):
    SH[i]=np.linalg.norm(F[i,:])    
SHrank=np.argsort(SH) # index of SH

# selected SH spanner
print "Selected SH spanners are:"
print SHrank[0:topk]   # Note: the index starts from 0, not from 1

## HAM_cluster_result
to_keep_index = np.sort(SHrank[topk:])
A_temp = A[to_keep_index, :]
A_temp = A_temp[:, to_keep_index]
HAM_labels_keep= np.asarray(ground_truth_labels)[to_keep_index]
allLabels = np.asarray(ground_truth_labels)

cluster_matrix = F
labelbook, distortion = kmeans(cluster_matrix[to_keep_index,:], c)
HAM_labels, dist = vq(cluster_matrix[to_keep_index,:], labelbook)

print "AMI"
print 'HAM: '+str(metrics.adjusted_mutual_info_score(HAM_labels, HAM_labels_keep.T[0]))

# classifify SHS using majority voting
predLabels = np.zeros(len(ground_truth_labels))
predLabels[to_keep_index] = HAM_labels+1
HAM_predLabels = label_by_neighbors(A,predLabels)
print 'HAM_all: '+str(metrics.adjusted_mutual_info_score(HAM_predLabels, allLabels.T[0]))

print "NMI"
print 'HAM: '+str(metrics.normalized_mutual_info_score(HAM_labels, HAM_labels_keep.T[0]))
print 'HAM_all: '+str(metrics.normalized_mutual_info_score(HAM_predLabels, allLabels.T[0]))

print "Entropy"
print 'HAM: '+str(avg_entropy(HAM_labels, HAM_labels_keep.T[0]))
print 'HAM_all: '+str(avg_entropy(HAM_predLabels, allLabels.T[0]))

# Save Results for SHII Metric Evaluation
fold = 1
datapath= os.getcwd() + '/SHII_metric/'
inputdata = 'datasets'

# Save selected top-k SHS
SHS = SHrank[0:topk] + 1
filename =  datapath+ inputdata +'/'+ adj_name + '_' + str(fold) + '_topk.txt'
shs_str = ' '.join([str(node) for node in SHS])
fout = open(filename,'w')
fout.write(shs_str+'\n')
fout.close()

# Save Edge Links (adjacency matrix)
filename =  datapath+ inputdata +'/'+  adj_name + '_' + str(fold) + '_edge_dir.txt'
sio.mmwrite(filename,A.astype('int'))
# Save Ground Truth Labels
filename =  datapath+ inputdata +'/'+  adj_name + '_' + str(fold) + '_label.txt'
save_labels(filename, ground_truth_labels)
