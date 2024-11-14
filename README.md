# 1 - Introduction

The code stuies the super-cluster phase transition in a random graph with all-to-all connections. The main calculation is provided by a C script, which exploits a tree-search algorithm to retrieve and store information about the individual clusters. A Python (Jupyter) notebook is provided which can either execute the C script or retrieve and analyzed previously stored data.

## 1.1 - The problem

We focus on a random graph of $N$ nodes, which can be pairwise connected by a link with a random (uniformly sampled) probability. Each node can be linked with up to $L_{\rm node} = N-1$ nodes, leading to a maximum amount of posible links given by $L_{\rm max} = N(N-1)/2$, where the factor $1/2$ avoids the double counting of identical links. In our implementation, we consider a Erdős–Rényi description of the system, which is defined from a number $N$ of nodes and a total number $L$ of distinct links, randomly (uniformly) distributed among all possible pairs. This treatment is closely related to the alternative Gilbert model, which defines the graph through a number $N$ of nodes and a probability $p_n$ for each node to be connected.

## 1.1.1 - Super-cluster phase transition
Given a number $N$ of nodes and a number $L$ of distinct, random links, we can define a cluster as an ensemble of nodes that are all connected by crossing one or more links (e.g. see the colored scheme in Fig. $a_1$). We can thus define the size $S$ of each cluster as the number of nodes contained in the cluster, and $S_{\rm max}$ as the size of the largest cluster. We define as super-cluster, one whose size is at least of the same order of $\sim \mathcal O (N)$. Similarly we can define the average number of connections per node as $\langle C \rangle=2L/N$, which is nothing more than the number of distinct links counted twice (since each link increases the number of connections of two distinct nodes) and normalized by the number of nodes.

In the thermodinamic limit of $N\to \infty$ it is easy to see that a phase transition is expected to occur when studying the maximal cluster size $S_{\rm max}$ as a function of the average number of connections $\langle C \rangle$. To show this, we notice that the probability that a node has exactly $C$ connections is given by the Bernoulli distribution $P(C) =  p^C (1-p)^{N-1-C} (N-1)!/ [C!(N-1-C)!]$, where we have the single-node connection probability $p = L/L_{\rm max} = 2L/[N(N-1)]$. In the limit of $N\to \infty$, the Bernoulli distribution simplifies to the Poisson distribution $P(C)\approx \{\langle L\rangle\}^C e^\{-\langle L\rangle\} / C!$, allowing to estimate the probability $p_{\rm in}$ that a (randomly chosen) node belongs to the largest cluster via the relation $1-p_\{\rm in\} = \sum_C P(C) \left(1-p_\{\rm in\}\right)^C = e^\{\langle C\rangle p_\{\rm in\}\}$. One than easily observe that if $\langle C\rangle<1$ the only valid solution is $p_\{\rm in\}=0$, meaning that in the limit of $N\to \infty$ the cluster size is $S_{\rm max}\sim p_\{\rm in\}N \sim \mathcal\{O\}(1)$, so no super-cluster exists. On the contrary, when $\langle C\rangle>1$, then a solution $p_\{\rm in\}\sim \mathcal O(1)$ exists, which is associated to a super cluster $S_{\rm max}\sim p_\{\rm in\}N \sim \mathcal\{O\}(N)$.



![random_graph_clusters](https://github.com/user-attachments/assets/d94fa873-b5e0-4257-9151-e8485bea0b25)

