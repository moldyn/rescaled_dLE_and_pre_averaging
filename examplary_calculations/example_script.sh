#!/bin/bash

name="LE_trajectory"

##### The following line creates a trajectory of 2*10^6 points based on the Markovian Langevin equation (LE) using a timestep of 0.02 fs, a temperature of 300 K and the files free_energy, mass, gamma, start as input for free energy, mass, friction and starting points

build/LE -start start -free free_energy -gamma gamma -mass mass -o ${name} -t 0.02 -T 300 -L 2000000 -d 1 -n 200

##### The following line uses the trajectory produced above and adds a column consisting of 1 indicating that the data does not consists of seperat short trajectories. The end of concatenated short trajectores would be indicated by a 0 in the added column

awk '!/^#/ {print $1,1}' ${name} > ${name}_inputdLE

#### The next line constructs a data-driven Langevin (dLE) trajectory based on our file ${name}_inputdLE

../dLE_normal/build/dLE -i ${name}_inputdLE -o ${name}.dLE.k300 -k 300 -L 2000000 -d 1

####  In the following, we produce a trajectory based on the generalized Langevin equation. Please note that the timestep of 0.002 ps needs to be so small to ensure convergence of the colored noise with decay time of 0.2 ps!

name="GLE_trajectory"

../Generalized_LE_const_fields/build/GLE -start start -free free_energy -gamma gamma -mass mass -tau tau -o ${name} -t 0.002 -s 10 -T 300 -L 2000000 -d 1 -n 200

##### Again, we add a column indicating that the produced trajectory represents a single simulation run...

awk '!/^#/ {print $1,1}' ${name} > ${name}_inputdLE

##### ... and apply the dLE to this file to get a new dLE-trajectory

../dLE_normal/build/dLE -i ${name}_inputdLE -o ${name}.dLE.k300 -k 300 -L 2000000 -d 1

### Analysing both dLE-trajectories will show that they differ quite drasticly with respect to the long time predictions although the input files didn't! For further reference please consider B. Lickert and G. Stock, "Modeling non-Markovian data using Markov state and Langevin models" mentioned in the READMEs of the different dLE-variations.
