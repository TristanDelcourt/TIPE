# TIPE: Factorisation implementations

## Presentation

This repo gives possible implementations of algorithms that are used to factor a large number into its prime factors. I have implemented 2 of those algorithms so far, these being:

* Dixon's algorithm
* The quadritc sieve (qsieve)
* The multipolynomial quadratic sieve (mpqs, as well as a parallelized version pmpqs)

You are able to use these algorithms from a single executable which I will present later

## General layout

The main parts of the algorithms implemented, as well as how to use them, are as follows:

* [dixon/](./c/dixon)
    * [dixon.c](./c/dixon/dixon.c)
* [qsieve/](./c/qsieve)
    * [qsieve.c](./c/qsieve/qsieve.c)
* [mpqs/](./c/mpqs)
    * [mpqs.c](./c/mpqs/mpqs.c)
    * [parallel_mpqs.c](./c/mpqs/parallel_mpqs.c)
* [Makefile](./c/Makefile)
* [main.c](./c/main.c)

In each algorithms' folder you can find the main part of how to find B-Smooth (defined in [[1](#sources)], page 5) relations according to that algorithm.

## Use of the algorithms

To build the program, run the [Makefile](./c/Makefile). After doing this, the factor executable will be created. Here are its arguments as well as what they do:

`--number %int%` (or `-n %int%`) 
This is the number to factor, the only required argument

`--type %arg%` (or `-t %arg%`) 
This sets what algorithm to use, currently these are 'dixon', 'qsieve', 'mpqs' or 'pmpqs' for parallelized mpqs

`--bound %int%` (or `-b %int%`) 
This sets the upper bound for the largest prime number to be used in the prime base

`--sieving_interval %int%` (or `-s %int%`) 
This sets the interval used for the sieveing of the quadratic sieve (will be ignored if dixon is used)

`--extra %int%` (or `-e %int%`) 
This sets the extra number of B-Smooth relations to find, by default this is 1

`-o %file_path%` 
This is the output file to send the results to if specified, otherwise it will print to the terminal

`--quiet` (or `-q`) 
Set this command to diseable all prints (useful for testing the algorithm without bloating the terminal)

## Sources

[[1](https://dspace.cvut.cz/bitstream/handle/10467/94585/F8-DP-2021-Vladyka-Ondrej-DP_Vladyka_Ondrej_2021.pdf?sequence=-1&isAllowed=y)] Quadratic sieve factorization algorithm, Bc. Ondˇrej Vladyka