# 1 - Introduction

The code stuies the super-cluster phase transition in a random graph with all-to-all connections. The main calculation is provided by a C script, which exploits a tree-search algorithm to retrieve and store information about the individual clusters. A Python (Jupyter) notebook is provided which can either execute the C script or retrieve and analyzed previously stored data.


## 1.1 - The problem

We focus on a random graph of $N$ nodes, which can be pairwise connected by a link with a random (uniformly sampled) probability. Each node can be linked with up to $L_{\rm node} = N-1$ nodes, leading to a maximum amount of posible links given by $L_{\rm max} = N(N-1)/2$, where the factor $1/2$ avoids the double counting of identical links. In our implementation, we consider a Erdős–Rényi description of the system, which is defined from a number $N$ of nodes and a total number $L$ of distinct links, randomly (uniformly) distributed among all possible pairs. This treatment is closely related to the alternative Gilbert model, which defines the graph through a number $N$ of nodes and a fixed probability for each node to be connected.


## 1.1.1 - Super-cluster phase transition

Given a number $N$ of nodes and a number $L$ of distinct, random links, we can define a cluster as an ensemble of nodes that are all connected by crossing one or more links (e.g. see the colored scheme in Fig. $a_1$). We can thus define the size $S$ of each cluster as the number of nodes contained in the cluster, and $S_{\rm max}$ as the size of the largest cluster. We define as super-cluster, one whose size is at least of the same order of $\sim \mathcal O (N)$. Similarly we can define the average number of connections per node as $\langle C \rangle=2L/N$, which is nothing more than the number of distinct links counted twice (since each link increases the number of connections of two distinct nodes) and normalized by the number of nodes.

In the thermodinamic limit of $N\to \infty$, it is easy to see that a phase transition is expected to occur when studying the maximal cluster size $S_{\rm max}$ as a function of the average number of connections $\langle C \rangle$. To show this, we notice that the probability that a node has exactly $C$ connections is given by the Bernoulli distribution $P(C) =  p^C (1-p)^{N-1-C} (N-1)!/ \[C!(N-1-C)!\]$, where we have the single-node connection probability $p = L/L_{\rm max} = 2L/\[N(N-1)\]$. In the limit of $N\to \infty$, the Bernoulli distribution simplifies to the Poisson distribution $P(C)\approx \{\langle L\rangle\}^C e^\{-\langle L\rangle\} / C!$, allowing to estimate the probability $p_{\rm in}$ that a (randomly chosen) node belongs to the largest cluster via the relation $1-p_\{\rm in\} = \sum_C P(C)  (1-p_\{\rm in\} )^C = e^\{-\langle C\rangle p_\{\rm in\}\}$. One than easily observe that if $\langle C\rangle<1$ the only valid solution is $p_\{\rm in\}=0$, meaning that in the limit of $N\to \infty$ the cluster size is $S_{\rm max}/N\sim p_\{\rm in\} \sim 0$, so no super-cluster exists. On the contrary, when $\langle C\rangle>1$, then a solution $p_\{\rm in\}\sim \mathcal O(1)$ exists, which is associated to a super cluster $S_{\rm max}/N\sim p_\{\rm in\} \sim \mathcal\{O\}(1)$.

In simple words, in a large $N\gg 1$ graph where each pair of nodes can be randomly connected, a super-cluster (whose size is comparable to $N$) appears if the nodes have, on average, at least one connection.

<br/>

![random_graph_clusters](https://github.com/user-attachments/assets/d94fa873-b5e0-4257-9151-e8485bea0b25)

<br/>

# 2 - The code

The code computes the size of the largest custer as a function of an increasing number of random links $L$, for a fixed number of nodes $N$. To do so, at each step two random nodes are selected and, if they are not already linked, a new link is added, which merges their respective clusters into a new larger cluster. The main algorithmic issue concerns the storage and retrieval of the information about the size and nodes of each cluster, as well as the update of this information when two clusters are merged.


## 2.1 - Simplest algorithm

The simplest treatment of the problem is described in Fig.s $a_1$ and $a_2$. Naively, one can indeed think of storing in each node two pieces of information: the nodes that are directly linked by a single link (which we will call neighbors), and a label identifying its cluster. Separately, a dictionary can be stored which associates each label with the number of nodes in that cluster. This treatment is exemplified in Fig. $a_1$, where each color identifies a different cluster. When two nodes are connected by adding a new link (as highlighted by the orange circle n Fig. $a_2$), it can happen that two clusters (blue and red) are then merged into a new cluster (purple). With this approach, merging two clusters requires updating the information in all the nodes, with a number of stepwise operations that can be as large as $S_\{\rm max\}$. Since this is computed up to $C=N$ steps, the simulation can roughly overheat up to $\sim N^2$, making it hard to simulate numbers of nodes as large as $N\sim 10^6\div10^7$.



## 2.2 - Tree-search algorithm

Instead of storing the global information in each single node of the cluster, this can be saved only into one single element, considered the head of the cluster (Fig.s $b_1$ and $b_2$). Then, each node conserves knowledge of only its neighbors (i.e. the node connected via a direct link, as shown by the solid black lines), and of which node introduced it into the cluster the first time (dashed, black lines). Then, when a new link between two nodes (say, node $A_i$ and node $B_j$) implies that two clusters must be merged (Fig. $b_2$), the number of required operations becomes much lower than the previous scenario. The process, indeed, only requires a tree search of the two heads $A_0$, $B_0$ of the clusters (red arrows in Fig. $b_2$), which only involves a small subset of the nodes. Then, a link is added between $A_i$ and $B_j$ (black, solid line), while the clusters are merged by making the head of the largest cluster (let's say $A_0$) the new head of the joint cluster, and then making the other head $B_0$ reference $A_0$ as the node that introduced it to the cluster, to keep the tree-like structure.

# 3 - Numerical results

In the plot below (Fig.s c and d), we report the simulation results. First, in Fig. c, we show the phase transition of the (normalized) maximal cluster size $S_{\rm max}/N$ as a function of the average node degree $\langle C\rangle = 2L/N$, which quantifies the average number of direct connections (i.e. neighbors) of each node. Each colored line identifies a scenario with an increasing number of nodes $N=10^2\div 10^7$. The plot exhibits a clear phase transition in the thermodinamic limit $N\to \infty$. Finally, in Fig. d, we study the finite-size effects, by plotting the value of $S_{\rm max}/N$ computed at the transition point $\langle C\rangle=1$, given different values of $N$. By fitting the data in a log-log scale, one can infer the scaling $S_{\rm max}(\langle C\rangle=1)/N\approx 2/N$, i.e. $S_{\rm max}(\langle C\rangle=1)\approx 2$.  


![cluster_percolation_plots](https://github.com/user-attachments/assets/16a70883-8553-4f57-8030-d9f2b06a5f52)
