GSAFQMC
=======

Example (simple) ground state AFQMC code for the Hubbard with some ED codes as basic checks. Will work code into a ground state chemistry code for edification. 

This repository consists of a ground-state CPMC, a ground-state ED, and a finite-temperature ED code to be used for the purposes of learning ground-state CPMC techniques. These codes are far from optimized - lapack should be used for matrix operations and the construction of states in the ED codes is lackluster at the moment. Nevertheless, these should be useful for understanding the basic algorithm. Future work will develop the CPMC code into useful for chemistry. 

The CPMC code takes in inputs in the afqmc.par file (the meanings of these inputs are at the bottom of the file) and outputs energies and the used parameters in the afqmc-parameters.dat and energy.dat files. 

The ED codes take in inputs from the ed.par files and return output energies in the energy.dat file. 

All codes may be compiled by typing make. The Makefile should be altered to suit your needs. 
