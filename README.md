# On the Probability of LU Factorization of Discrete Random Matrices

This repository contains the code and computational results for the project, "On the Probability of LU Factorization of Discrete Random Matrices."

## Authors and Affiliation

  * Samuel Orellana Mateo (Duke University)
  * John Urschel (Massachusetts Institute of Technology)
  * Nicholas West (Massachusetts Institute of Technology)

This work was conducted as part of the MIT Summer Research Program (MSRP) in collaboration with the MIT Office of Graduate Education and the MIT Department of Mathematics.

## Usage

First, ensure the packages `StaticArrays` and `Glob` are installed in your Julia environment.

To run the high-performance generation script (uses multithreading and optimizations):
```
julia --threads=auto -O3 --check-bounds=no main.jl
```
To convert the generated binary data into readable JSON format:
```
julia readable.jl
```
To compute the upper bound on the probabilities:
```
julia upper_bound.jl
```

## Acknowledgments and License

This project utilizes algorithms and logic from the nauty and Traces package (v2.9.3).
Original software by Brendan McKay and Adolfo Piperno.
Licensed under Apache License 2.0.