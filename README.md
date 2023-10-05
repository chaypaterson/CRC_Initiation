hopscotch
=========

A Monte Carlo simulation suite for [our 2020 paper on initiation of colorectal
adenocarcinoma](https://www.pnas.org/doi/10.1073/pnas.2003771117).

Note that the indices for mutation and LOH are switched relative to the main
paper. For both APC and TP53, the indices of the different states in
params.h and hopscotch.cpp should be read as follows:

    0 = wild-type (tumour suppressor +/+)
    1 = mutated (t.s. +/-)
    2 = loss of chromosome/LOH (t.s. +/-)
    3 = one copy mutated, other chromosome lost (t.s. -/-)
    4 = both copies mutated (t.s. -/-).

1. How to build
---------------

Dependencies
------------

The code is C++11 compliant, and depends on the GNU Scientific Library for
some of the random number generation, libpthread for parallelisation, and 
math.h for other functions. 

On Debian/Ubuntu, the math.h and libpthreads libraries should be pre-installed.

The following gsl libraries can be installed with: 

apt install libgsl-dev libgsl23 libgslcblas0

They are present in the default repositories on Debian 8+ and Ubuntu 18+.

Compiling
---------

On Debian and Ubuntu, the code should then compile with

g++ hopscotch.cpp --std=c++11 -lgsl -lgslcblas -lm -o OUTPUTEXE -pthread

On MacOS, it should compile with

g++ hopscotch.cpp --std=c++11 -pg -lgsl -o OUTPUTEXE 

Note: on MacOS, you will have to add additional arguments to point the
linker to your local installation of pthreads and GSL. This will depend on
the library directory in your Homebrew or Fink installation.

2. How to run
-------------

Required arguments:

At a minimum, the program needs an output file. The default mode, if no
mode is specified, is tau-leaping.

At most one of the following modes should be provided:
--default : tau-leaping
--gillespie : Gillespie SSA
--lineages : Lineage-tracking tau-leaping

If more than one mode is specified, the program will exit with an error.

The output file is given with the argument:

-o OUTPUTFILE : OUTPUTFILE is the string that will be the output CSV file
                name.

If no output file is specified, the program will also exit with an error.

Other optional run-time options:

--cores N: N = number of threads
--withkey : prints the key in the first line of the output CSV file
--polite : with parallelisation, sets the number of threads to 50% of
           available cores
--host N : N = seed for RNG (can be integer e.g. output of date +%s, or string 
            e.g. hostname)

Any other options parsed in hopscotch.h are defunct, and should be ignored.

3. Example
----------

The most detailed lineage-tracking data in the paper was generated with

time nice ionice ./hopscotch --lineages --host 0 --cores 80 -o noFAP-1e4x80runs-bapc0.20-bkras0.07-bboth0.27-lineagetauleaping-tau0.01-seed0-newmultipliers-skipping.csv

Note: that tau has to be set inside hopscotch.cpp. You will need to change
line 736:

double dt = 0.01;

to set tau to some other value than 0.01.
