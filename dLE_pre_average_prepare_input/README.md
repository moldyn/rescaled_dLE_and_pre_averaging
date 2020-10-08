Script used to "pre-average" input data to accelarate data-driven Langevin
calculations. For further information see

*B. Lickert and G. Stock, "Modeling non-Markovian data using Markov state and 
Langevin models", in preparation*

For practical information please compile the program by running 

```
g++ prepareinput.cpp -o prepareinput -Wall -pedantic
```

in the shell.

The executable program will be stored in the folder build. By running

```
build/prepareinput -h
```

the program displays the available input options in detail together with some 
additional information. 

Intepretation of results: Please consult the reference above for
information on the interpretation of the differnt columns of the output file.
