# UTB Symbolic Python

(This file last edited 06-01-2019)

This repository is a space for fleshing out ideas for code that will be used in a quantum control package meant to perform the types of optimizations in the repository Robust-quantum-control. The goal is to build symbolic tools using SymPy that allow a user to specify their control problem symbolically (i.e. structure of the DE, the Dyson terms they wish to compute), and result in automatically generated code that computes these terms, and their derivatives, in a way that minimizes the number of matrix multiplications required (in, e.g. exponentiation). 

By it's nature, this repository is not meant to be usable, but is just a convenient place to figure out what will work, while keeping things sharable. 