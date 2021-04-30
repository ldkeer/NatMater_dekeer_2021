This GitHub repository contains additional input/output files accompanying the manuscript

**Computational prediction of the molecular configuration of three-dimensional network polymers**

by Lies De Keer, Karsu I. Kilic, Paul H.M. Van Steenberge, Lode Daelemans, Daniel Kodura, Hendrik Frisch, Karen De Clerck, Marie-Françoise Reyniers, Christopher Barner-Kowollik, Reinhold H. Dauskardt, and Dagmar R. D’hooge.

The authors declare that the complete algorithm of the modeling framework is fully available in the Supplementary Information of this contribution. The authors specifically refer to all implementation steps included in Supplementary Figure 1-4 and their discussion in the Supplementary Discussion complemented with information already available in the open literature.

Additional information concerning the data is available upon request from the authors. Please send a mail to Dagmar.Dhooge@UGent.be for more information.

**Table of Contents**

_Gillespie algorithm:_    Stochastic Simulation Algorithm (SSA) to perform kinetic Monte Carlo simulations; contains Fortran file 
                          (MonteCarloAlgorithm.for) and illustration of algorithm (GillespieAlgorithm.bmp)
                          
_LAMMPS_input:_     Datafiles include the information related to the initial computation box dimensions, atom types and masses, together with bond, angle and dihedral potential parameters; the initial topology information including the x, y, z coordinates of the atoms together with the bond pairs, angle triplets and dihedral quadruplets. There are three datafiles (data.SiO2_2000_1.1, data.SiO2_2000_1.2 and data.SiO2_2000_1.3) and each corresponds to a different initial condition with different sets of atom coordinates, bond pairs, angle triplets and dihedral quadruplets. The main input file (in.SiO2) has the information about the simulation dynamics which includes the types of bonded and non-bonded interactions used, simulated annealing parameters such as the time step, number of time steps, temperature and pressure, and the type of ensembles used in the different sections of the simulation. First a soft potential is applied on the initial configuration to prevent possible overlaps. Afterwards, simulated annealing is applied with Stillinger-Weber potentials to generate the network, which has two different steps in which the system temperature and pressure are decreased with NPT dynamics. OCS.sw includes empirical Stillinger-Weber potential parameters for each Si-O triplet. This file is called inside the input file during the simulated annealing part of the simulation.

_LAMMPS_output:_    Output files (SiO2_160000_1.1.xyz, SiO2_160000_1.2.xyz and SiO2_160000_1.3.xyz) include the final atom coordinates of the equilibrated silica networks. Each row in the output files includes the atom id, atom type and x, y, z coordinates of the corresponding atom.
                          
_MTS:_    Algorithm to calculate minimal spanning tree of graphs out of topology matrix; contains Fortran file (Main.for), input file (rede.dat) and output file (ARVORE.res)

_3DVisualization.for:_    Fortran file to identify network molecules using depth-first strategy with molecular dynamics simulation results as input; to create input file for Gaussview

_BinaryTreeOperations.txt:_     Pseudo-code related to operations in a binary tree (for more information, reader is refered to: P.H.M. Van Steenberge, D.R. D’hooge, M.-F. Reyniers, G.B. Marin, Improved kinetic Monte Carlo simulation of chemical composition-chain length distributions in polymerization processes, Chemical Engineering Science, 110, 2014, 185-199)

_Gaussview.dat:_    Example input file for Gaussview for 3D visualization of a network molecule
                          
_Gephi.csv:_    Example input file for Gephi for 2D visualization of network molecules out of topology matrix

_dijkstra_openmp.for:_    Fortran subroutine of Dijkstra graph distance algorithm using openMP application program interface for parallelization purposes to determine shortest path between two functional groups in a network molecule, necessary to apply the distance rule

_searchCycles_list.m:_    Matlab script to calculate molecular pore size distribution out of topology matrix
