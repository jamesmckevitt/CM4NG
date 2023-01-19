# Computing the gravitational potential on nested meshes using the convolution method
Convolution method for nested grids (CM4NG)

This code computes the gravitational potential on nested meshes 
using the convolution theorem. The method is thoroughly described 
in Vorobyov, McKevitt, et al., A&A, 2023, in press.

## Copyright statement

Computing the gravitational potential on nested meshes using the 
convolution method Copyright © 2023, Eduard Vorobyov and James 
McKevitt. All rights reserved. Anybody interested in using this 
code should contact the authors.

- eduard.vorobiev@univie.ac.at
- james.mckevitt@univie.ac.at

## Usage instructions

A Makefile is used, which for the GPU version requires some 
configuration for the specific destination GPU and location of 
intel compiler libraries.

The compilation is done by typing 'make'. The programme can be run
using the command 'make run'. Cleaning the directory is done with 
'make clean'.

Use 'mesh.f' to set: 
1) the number of grid zones per dimension 'N',
2) the number of nested grids 'Ndepth',
3) the desired density configuration (four choices are available: 
homogeneous sphere, homogeneous oblate ellipsoid, wide-separation 
binary, and non-homogeneous sphere with a density structure after 
Wand & Yen (2020)),
4) the number of OpenMP threads 'Nthreads'.

Use 'main.f' to choose the computation method: basic or advanced. 
The advanced method correctly computes the dipole moments of a 
close binary across each pair of neighboring nested grids (see 
Sect. 7).

In the 'main.f' file one can also choose either to use the 
convolution theorem or direct summation when computing the input 
from the doubled mesh in the advanced method (see the Appendix).

## Additional notes

The CPU code has been verified with Intel 19.1.3.304. The GPU 
version has been verified with Intel 19.1.3.304 and NVIDIA HPC SDK 
22.5-gcc-11.2.0-h6et62m, and has been tested on GPUs with CUDA 
version 11.6 and compute capabilities 8.0 and 8.6. Details of the 
specific hardware used for performance tests is in Sect. 6.

An external library 'libacml_mp.a' (from the ACML library package) 
is used to run the fft subroutines on CPUs. This is contained in 
this distribution.
