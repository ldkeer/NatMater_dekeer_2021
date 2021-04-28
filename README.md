This GitHub repository contains additional input/output files accompanying the manuscript

**Computational prediction of the molecular configuration of three-dimensional network polymers**

by Lies De Keer, Karsu I. Kilic, Paul H.M. Van Steenberge, Lode Daelemans, Daniel Kodura, Hendrik Frisch, Karen De Clerck, Marie-Françoise Reyniers, Christopher Barner-Kowollik, Reinhold H. Dauskardt, and Dagmar R. D’hooge.

The authors declare that the complete algorithm of the modeling framework is fully available in the Supplementary Information of this contribution. The authors specifically refer to all implementation steps included in Supplementary Figure 1-4 and their discussion in the Supplementary Discussion complemented with information already available in the open literature.

Additional information concerning the data is available upon request from the authors. Please send a mail to Dagmar.Dhooge@UGent.be for more information.

**Table of Contents**

_Gillespie algorithm:_    Stochastic Simulation Algorithm (SSA) to perform kinetic Monte Carlo simulations; contains Fortran file 
                          (MonteCarloAlgorithm.for) and illustration of algorithm (GillespieAlgorithm.bmp)
                          
_MTS:_  Algorithm to calculate minimal spanning tree of graphs out of topology matrix; contains Fortran file (Main.for), input file (rede.dat) and output file (ARVORE.res)

_3DVisualization.for:_  Fortran file to identify network molecules using depth-first strategy with molecular dynamics simulation results as input; to create input file for Gaussview

_BinaryTreeOperations.txt:_  Pseudo-code related to operations in a binary tree (for more information, reader is refered to: P.H.M. Van Steenberge, D.R. D’hooge, M.-F. Reyniers, G.B. Marin, Improved kinetic Monte Carlo simulation of chemical composition-chain length distributions in polymerization processes, Chemical Engineering Science, 110, 2014, 185-199)

_Gaussview.dat:_  Example input file for Gaussview for 3D visualization of a network molecule
                          
_Gephi.csv:_  Example input file for Gephi for 2D visualization of network molecules out of topology matrix

_dijkstra_openmp.for:_  Fortran subroutine of Dijkstra graph distance algorithm using openMP application program interface for parallelization purposes to determine shortest path between two functional groups in a network molecule, necessary to apply the distance rule

_searchCycles_list.m:_  Matlab script to calculate molecular pore size distribution out of topology matrix
