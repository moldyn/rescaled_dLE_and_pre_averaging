Script to propagate the generalized Langevin equation based on a exponential memory kernel.
The memory kernel is assumed to follow K_{i,j}=\Gamma_{i,j}/\tau_i * e^{-t/tau_i} with tau_{i}
being the decay time of coordinate i.

Constant friction matrix \Gamma and constant mass are assumed. For further information
please compile the program using CMake by running 

```
cmake -H. -Bbuild
cmake --build build -- -j3
```

in the shell. **The Eigen package is needed for successful compilation.**

The executable program will be stored in the folder build. By running

```
build/GLE -h
```

the program displays the available input options in detail together with some information
on the used integrator. 

Unit convention: The script assumes that the free energy is given in units of k_BT. In addition,
k_BT is measured in ps^-1, 300 K refer to k_BT=38 ps^-1. Masses are assumed to be given in ps.
A mass of 26u corresponds to 400 ps. Friction \Gamma is assumed to be without any unit, the noise 
amplitude K=\sqrt(2k_BT\Gamma) follows from \Gamma and k_BT.

**Providing the different fields in wrong units leads to wrong results!!!**
