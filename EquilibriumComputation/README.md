This repository contains the files necessary for the computation of the equilibrium. The file main.cpp contains all the necessary function in c++. The makefile assumes that we can use a given number of cores to parallelize. In order to run the program we need to have boost and gsl library. Nlopt should be installed with the commands:

```
make nloptinst
```

In order to execute silently:

```
make executesilent
```

The output will be a file called nohupEXECUTE.out containint the output from the terminal and a csv file called NashEquilibrium containint the relevant strategies and payoffs in NEQ. In order to run it and see the output:

```
make execute
```
