# SU-2-Wilson-Action
In this project I'm trying to simulate an SU(2) gauge group on a 3+1 Euclidean lattice. Specifically, I'm trying to recreate one of the graphs in Creutz's book, "Quarks, Gluons, and Lattices," for confirmation that my simulation is correct. From there, I want to find the eigenvalues of the Wilson-Dirac operator (generated from the lattice configuration) using the Hermitian Lanczos algorithm and then a reguiar eigenvalue finding algorithm. In this project I use the GSL to do most of the complex number/matrix operations.  

    Installation instructions: make sure that the GNU scientific library (GSL) is properly installed before running the make file. 
