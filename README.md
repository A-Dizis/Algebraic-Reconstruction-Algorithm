# Algebraic-Reconstruction-Algorithm
This is an application of ART algorithm in 2-D using multithreading and a new approach on the calculation of the projection matrices.

For the user to get a grasp of the internal workings of this algorithm, it is highly recommended to, initially skip through the paper "A tutorial on ART (Algebraic Reconstruction Technique)" by R. Gordon.

In this paper the projection matrices are not calculated with accuracy, as they are given bit values 0 or 1 in each position. Here we calculate the values with great accuracy, using lattice points for the estimation of the real area of each position. 
