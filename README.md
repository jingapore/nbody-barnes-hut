# 1. TLDR
I present findings from a parallel implementation of Barnes-Hut n-body simulation, where we achieve up to 6x speedup from 1 process to 8 process.

# 2. Design behind Parallel Implementation
## 2.1 Central Problem of Parallelisng Barnes-Hut
There are two stages to Barnes-Hut: (a) constructing a tree, and (b) computing forces on each particle from the tree. Parallelising force computation in stage (b) is trivial. However, it is not clear how we can parallelise the construction of the tree in stage (a).

Assuming stages (a) and (b) take roughly equal time (given that their theoretical time complexity is O(n log n) each), Amdahl’s law indicates that going from 1 processor to 8 processor can only result in 1.78x speedup (with only 50% of the entire programme parallelised).

However, this lab has achieved at 6x speedup when going from 1 process to 8 processors. We achieved this using a parallel setup defined in Salmon (1991), which this report will explore in some depth.

The key idea from Salmon (1991) is to construct locally sufficient subtrees within each processor. Each processor will be responsible for particles within a spatial segment, and will construct their respective subtrees. However, each subtree will not yet be sufficient for force computation. The subtree is only locally sufficient after a round of data exchange between processors, facilitated by MPI as the message passing model. 

## 2.2 Domain Decomposition
For each iteration of a Barnes-Hut time-step, we ask ourselves: how can we decompose the domain across different processors, such that each processor is assigned fairly equal work? If we fail to do so, load imbalances will reduce gains from parallelism; some processors will be idle while other processors are performing effortful force computations.

The naive approach is to decompose the space into equal segments of n_processors, and allocate each processor to a segment. This would have resulted in a uniform grid for 2D, or a cube for 3D.

However, this naive approach encounters two problems. First, certain spatial segments may have greater concentration of particles. Second, even if each processor handles the same number of particles, the force computation workload might not be evenly spread. This unevenness happens because each update to a particle requires different traversal depths of the Barnes-Hut tree.

The solution proposed in Salmon (1991), is to use orthogonal recursive bisection that is sensitive to the workload that had just been executed by the particle in the previous iteration. This is intuitively sound. Assuming a small enough dt, the update to each particle should be incremental, and so the work required in the current iteration should not deviate materially from the work in the previous iteration. Figure 1 illustrates this idea.

<p align="center"><b>Figure 1: Comparison of Naive Bisection v.s. Work-sensitive Bisection</b></p>
<p align="center">
    <img width="66%" src="img/fig_1.png">
</p>

Each processor establishes log n number of partners, on the other side of a bisection. We index each partner using bitwise operations. The usefulness of each processor having a partner across each bisection will be apparent in sub-section 2.3. For now, we already have an algorithm for orthogonal recursive bisection, which enables even load balancing across processors.

Before we proceed to the next section, we should be aware that there are two blocking operations across processors.
- First, each processor handles a segment of bisection, for every iteration of a time step. However, each processor is unable to have a global view of the work required within each spatial segment. So our program requires message passing whenever it considers the suitability of a bisection split.
- Second, after every bisection, the processors assigned to each bisection may not have all the particles they are responsible for. This requires an exchange of bodies with partner processors, which are at the other side of the bisection.

# References
- Barkman, P. (2017) et al. N-body. GitHub. https://github.com/barkm/n-body
- Salmon, J. K. (1991). Parallel hierarchical N-body methods (Doctoral dissertation). California Institute of Technology.
- Salmon, J. K. & Warren, M (1993). A Parallel Hashed Oct-Tree N-Body Algorithm. California Institute of Technology & Los Almos National Laboratory.
