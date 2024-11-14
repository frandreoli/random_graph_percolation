# 1 - Introduction

The code stuies the super-cluster phase transition in a random graph with all-to-all connections. The main calculation is provided by a C script, which exploits a tree-search algorithm to retrieve and store information about the individual clusters. A Python (Jupyter) notebook is provided which can either execute the C script or retrieve and analyzed previously stored data.

## 1.1 - The problem

We focus on a random graph of $N$ nodes, which can be pairwise connected by a bridge with a random (uniformly sampled) probability. Each node can be bridged with up to $B_{\rm node} = N-1$ nodes, leading to a maximum amount of posible bridges given by $B_{\rm max} = N(N-1)/2$, where the factor $1/2$ avoids the double counting of identical bridges.
