Efficient Global Optimal Resource Allocation in Non-Orthogonal Interference Networks
==================

This code package is related to the following scientific article:

Bho Matthiesen and Eduard Jorswieck, "[Efficient Global Optimal Resource Allocation in Non-Orthogonal Interference Networks](https://arxiv.org/abs/1812.07253)," submitted to IEEE Transactions on Signal Processing.


## Abstract of Article

Many resource allocation tasks are challenging global (i.e., non-convex) optimization problems. The main issue is that the computational complexity of these problems grows exponentially in the number of variables instead of polynomially as for convex optimization problems. However, often the non-convexity stems only from a subset of variables. Conventional global optimization frameworks like monotonic optimization or DC programming [1] treat all variables as global variables and require complicated, problem specific decomposition approaches to exploit the convexity in some variables [2]. To overcome this challenge, we develop an easy-to-use algorithm that inherently differentiates between convex and non-convex variables, preserving the polynomial complexity in the number of convex variables. Another issue with these widely used frameworks is that they may suffer from severe numerical problems. We discuss this issue in detail and provide a clear motivating example.The solution to this problem is to replace the traditional approach of finding an ε-approximate solution by the novel concept of ε-essential feasibility. The underlying algorithmic approach is called successive incumbent transcending (SIT) algorithm and builds the foundation of our developed algorithm. A further highlight of our algorithm is that it inherently treats fractional objectives making the use of Dinkelbach's iterative algorithm obsolete. Numerical experiments show a speed-up of five orders of magnitude over state-of-the-art algorithms and almost four orders of magnitude of additional speed-up over Dinkelbach's algorithm for fractional programs.

References:

1. H. Tuy, Convex Analysis and Global Optimization, ser. Springer Optimization and Its Applications. Springer, 2016.
2. B. Matthiesen and E. A. Jorswieck, "Weighted sum rate maximization for non-regenerative multi-way relay channels with multi-user decoding," in Proc. IEEE 7th Int. Workshop Comput. Adv. Multi-Sensor Adaptive Process. (CAMSAP), Curaçao, Dutch Antilles, Dec. 2017.

## Requirements & Contents of the Code Package

This code was compiled and tested with GNU Make 3.82 and GCC 7.3. It requires [Gurobi](http://www.gurobi.com/) and [Mosek](https://www.mosek.com/) which are both available free of charge for reserach in academic institutions. The simulations for the article cited above were done with Gurobi 8.0.1 and Mosek 8.1.0.34.

The complete source code is available in `code/`. Prior to compilation, please update the variables in the Makefile according to your needs. The simulations were conducted on TU Dresden's HPC. Sample [SLURM](https://www.schedmd.com/) `sbatch` files are available in `slurm/`. All input data is stored in `data` and raw results plus evaluation scripts are in `results/`.

## Acknowledgements

This research was supported in part by the Deutsche Forschungsgemeinschaft (DFG) in the [Collaborative Research Center 912 "Highly Adaptive Energy-Efficient Computing."](https://tu-dresden.de/ing/forschung/sfb912) and under grant number JO 801/24-1.

We thank the Center for Information Services and High Performance Computing (ZIH) at TU Dresden for generous allocations of computer time.


## License and Referencing

This program is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

