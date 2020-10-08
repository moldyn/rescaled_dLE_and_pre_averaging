Script to apply the data-driven Langevin equation without modifications. For further
information see

*N. Schaudinnus, B. Lickert, M. Biswas, and G. Stock, "Global Langevin model of 
multidimensional biomolecular dynamics", J. Chem. Phys. 145, 184114 (2016)*

or

*N. Schaudinnus, B. Bastian, R. Hegger,and G. Stock, "Multidimensional Langevin 
modeling of nonoverdamped dynamics", Phys. Rev. Lett. 115, 050602 (2015)*

For practical information please compile the program using CMake by running 

```
cmake -H. -Bbuild
cmake --build build -- -j3
```

in the shell. **The Eigen package is needed for successful compilation.**

The executable program will be stored in the folder build. By running

```
build/dLE -h
```

the program displays the available input options in detail together with some 
additional information. 

Intepretation of results: Please consult the references above for
information on the interpretation of the differnt columns of the output file.
