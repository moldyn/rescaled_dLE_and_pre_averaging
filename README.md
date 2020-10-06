# Quick Start

It is highly recomended to read the references given below. Besides this, every subfolder contain its own README-file
explaining how to complile the individual programs.

Every program can be executed with

program -h 

to get information about the input options and some general aspects of the code.

Besides this, the folder 

examplary_calculations

contains everything to do the first dLE calculations. It should be enough to get the logic behind the different
programs so that their application should not be too hard.

# Description of different folders

The different folders contain programs to construct trajectories from Langevin models. First, there are
three folders for Langevin models where free energy, friction and system mass are already known

Markov_LE_const_fields: This folder contains a program which uses a given free energy, a given constant
friction and a given constant mass to produce a (possible high-dimensional) Markovian Langevin trajectory.
Markov_LE_vary_fields: Here, friction and mass used to produce the trajectory can show some dependence on
the coordinate x.
Generalized_LE_const_fields: Here, a trajectory can be produced based on the generalized Langevin equation.
The memory kernel is assumed to decay mono-exponential and it is assumed to be independent of coordinate x.
The mass is assumed to be independent of x as well.

Then, there are three folders considering the data-driven Langevin equation (dLE) which constructs the different
forces (free energy, friction and stochastic force) based on given input data. For reference please have a look at

N. Schaudinnus, B. Lickert, M. Biswas, and G. Stock, "Global Langevin model of 
multidimensional biomolecular dynamics", J. Chem. Phys. 145, 184114 (2016)

or

N. Schaudinnus, B. Bastian, R. Hegger,and G. Stock, "Multidimensional Langevin 
modeling of nonoverdamped dynamics", Phys. Rev. Lett. 115, 050602 (2015)

The three folders are

dLE_normal: Produces a new trajectory based on input data.
dLE_normal_reflect: Allows to reflect the Langevin dynamics during propagation at an upper and a lower border.
Each coordinate x_i can have different borders
dLE_normal_testmodel: Constructs the noise observed in a given input trajectory in the Markovian Langevin framework.
Can be used to verify the model.

The next three folders contain a variation of the dLE where the friction can be rescaled to modify the model
dynamics. For consistency reasons, the integration time step needs to be rescaled as well. Further information
can be found at

B. Lickert and G. Stock, "Modeling non-Markovian data using Markov state and 
Langevin models", in preparation

The three folders are

dLE_rescaled: Does conceptionally the same as dLE_normal allowing for the mentioned rescaling.
dLE_rescaled_reflect: Does conceptionally the same as dLE_normal_reflect allowing for the mentioned rescaling.
dLE_rescaled_testmodel: Does conceptionally the same as dLE_normal_testmodel allowing for the mentioned rescaling

Finally, there are four folders considering a dLE implementation which works with pre-averaged input data, i.e.,
this implementation runs faster than the other ones. Please have, again, a look at 

B. Lickert and G. Stock, "Modeling non-Markovian data using Markov state and 
Langevin models", in preparation

for reference.

The four folders are

dLE_pre_average_prepare_input: Contains a program which does the mentioned pre-averaging. 
dLE_pre_average: Does conceptionally the same as dLE_normal using preaveraged data.
dLE_pre_average_testmodel: Does conceptionally the same as dLE_normal_testmodel using preaveraged data.
dLE_pre_averaged_rescaled: Does the same as dLE_rescaled but uses pre-averaged data, i.e., it is possible 
to rescale the friction (and the integration time step).

The final folder examplary_calculations contains everything needed to do the first dLE runs (based on dLE_normal)
after producing input data with Markov_LE_const_fields and/or Generalized_LE_const_fields. It should help to understand
the logic behind the dLE framework.
