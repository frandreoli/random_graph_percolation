# 1 - Introduction

The code stuies the super-cluster phase transition in a random graph with all-to-all connections. The main calculation is provided by a C script, which exploits a tree-search algorithm to retrieve and store information about the individual clusters. A Python (Jupyter) notebook is provided which can either execute the C script or retrieve and analyzed previously stored data.

## 1.1 - The problem

We focus on a random graph of $N$ nodes, which can be pairwise connected by a bridge with a random (uniformly sampled) probability. Each node can be bridged with up to $B_{\rm node} = N-1$ nodes, leading to a maximum amount of posible bridges given by $B_{\rm max} = N(N-1)/2$, where the factor $1/2$ avoids the double counting of identical bridges. In our implementation, we consider a Erdős–Rényi description of the system, which is defined from a number $N$ of nodes and a total number $B$ of distinct bridges, randomly (uniformly) distributed among all possible pairs. This treatment is closely related to the alternative Gilbert model, which defines the graph through a number $N$ of nodes and a probability $p_n$ for each node to be connected.

## 1.1.1 - Super-cluster phase transition
Given a number $N$ of nodes and a number $B$ of distinct, random bridges, we can define a cluster as an ensemble of nodes that are all connected by crossing one or more bridges (e.g. see the colored scheme in Fig. $a_1$). ![random_graph_clusters](https://github.com/user-attachments/assets/b91d77c0-82da-47fa-b299-592749ed0c87)
