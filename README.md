# 1 - Introduction

The code stuies the super-cluster phase transition in a random graph with all-to-all connections. The main calculation is provided by a C script, which exploits a tree-search algorithm to retrieve and store information about the individual clusters. A Python (Jupyter) notebook is provided which can either execute the C script or retrieve and analyzed previously stored data.

## 1.1 - The problem

We focus on a random graph of $N$ nodes, which can be pairwise connected by a bridge with a random (uniformly sampled) probability. Each node can be bridged with up to $B_{\rm node} = N-1$ nodes, leading to a maximum amount of posible bridges given by $B_{\rm max} = N(N-1)/2$, where the factor $1/2$ avoids the double counting of identical bridges. In our implementation, we consider a Erdős–Rényi description of the system, which is defined from a number $N$ of nodes and a total number $B$ of distinct bridges, randomly (uniformly) distributed among all possible pairs. This treatment is closely related to the alternative Gilbert model, which defines the graph through a number $N$ of nodes and a probability $p_n$ for each node to be connected.

## 1.1.1 - Super-cluster phase transition
Given a number $N$ of nodes and a number $B$ of distinct, random bridges, we can define a cluster as an ensemble of nodes that are all connected by crossing one or more bridges (e.g. see the colored scheme in Fig. $a_1$). We can thus define the size $S$ of each cluster as the number of nodes contained in the cluster, and $S_{\rm max}$ as the size of the largest cluster. Similalry we can define the average number of connections per node as $<C>=2B/N$, which is nothing more than the number of distinct bridges counted twice (since each bridge increases the number of connections of two distinct nodes) and normalized by the number of nodes.

In the thermodinamic limit of $N\to \infty$ it is easy to see that a phase transition is expected to occur when studying the maximal cluster size $S_{\rm max}$ as a function of the avrage number of connections $<C>$. To show this, we notice that the probability that a node has exactly $C$ connections is given by the Bernoulli distribution $P(C) =  p^C (1-p)^{N-1-C} (N-1)!/ [C!(N-1-C)!]$, where we have the single-node connection probability $p = B/B_{\rm max} = 2B/[N(N-1)]$. In the limit of $N\to \infty$ the Bernoulli distribution simplifies to the Poisson distribution $P(C)\approx <B>^C e^{-<B>} / C!$




![random_graph_clusters](https://github.com/user-attachments/assets/d94fa873-b5e0-4257-9151-e8485bea0b25)

