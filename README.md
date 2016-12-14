# CS466 Final Project

## 1) Background
In order to apply the algorithm, we first created the benchmark by creating 70 data sets with different parameters for ICPC, ML, and SC. Each dataset has the motif, motiflength, sequences, and planted sites stored in separate files. We generated all these random sequences and motifs and planted the sites randomly as well with the different parameters, repeated up to 10 times with 7 different parameter combinations.
## 2) Our Algorithm
We are using a Greedy Motif Search. First we are finding the best alignment between 2 sequences by trying every possible alignment and maximizing the score. Next, we are incrementally adding another sequence and trying every possible alignment with the new sequence against our current alignments to again maximize the score. 

This results in a running time of O(L^2 + (t-2)L ) where t is the number of sequences and L is the total number of starting positions in each sequence.
