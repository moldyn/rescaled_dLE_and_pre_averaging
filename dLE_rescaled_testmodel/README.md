Script to calculate the noise trajectory of some given data in the framework of
the data-driven Langevin equation **with rescaling of friction and integration 
timestep.** For further information see

*B. Lickert and G. Stock, "Modeling non-Markovian data using Markov state and 
Langevin models", in preparation*

For practical information please compile the program using CMake by running 

```
cmake -H. -Bbuild
cmake --build build -- -j3
```

in the shell. **The Eigen package is needed for successful compilation.**

The executable program will be stored in the folder build. By running

```
build/dLE_rescaled_testmodel -h
```

the program displays the available input options in detail together with some 
additional information. 

Intepretation of results: Please consult the reference above for
information on the interpretation of the differnt columns of the output file.
