# NetworkTopology

## Requirements

1. MATLAB (versions 2018b and 2019b were used in creating codes)
2. MATCONT (version 6.11)

## Getting Started

Clone the repository:
    
    git clone https://github.com/Shvartsman-Lab/NetworkTopology.git
    
## Running Codes

neywork_oscillations.m is the main script for this project.  After setting the value of the parameter delta to be studied, as well as the network topology type and number of cells, the system of ODEs are run.  The outputs are stored as the X and Y values for each nurse cell under consideration, and three different plots are shown: the X,Y phase plane for each nurse cell, the time evolution of X_i for each nurse cell, and a 3D phase plot of X_2, X_3, X_4 (much of the analysis in the paper is done on networks with three nurse cells, but this can be edited as necessary).

unp.mat and unpy.mat are datasets containing points on the (X,Y) limit cycle in the absence of inhibitor.

simple_cyc_ode.m is a function that defines the system of equations for the evolution of X_i and Y_i on whatever prescribed network is chosen for analysis.

build_adj.m is a function that builds the tree, line, or star matrix based on the number of cells desired (for the tree case, only 2,4,8, or 16 cell networks can be generated).
