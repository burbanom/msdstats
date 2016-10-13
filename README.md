# msdstats
Calculate statistics of mean square displacements and diffusion coefficients from molecular dynamics simulations.

The package consists of a python script and a module written in fortran.
To generate the fortran module it is necessary to compile it as follows:

f2py -c -c calcmsds msdconf.f90

The calcmsds.so should be placed in a folder that is defined in your $PATH/$PYTHONPATH variable.
