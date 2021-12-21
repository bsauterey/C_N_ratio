# C_N_ratio

To run for compilation: "gcc -o code Shared_Code_lat_input.c -fopenmp"
After compilation, this version of the code compute the eco-evolutionary dynamics
of the plankton community in one location of the North-Atlantic transect, and takes
the latitude of that location as argument (i.e., type "./code lat"). It can be
modified to do the same for multiple locations (i.e., by setting Spaceres > 0).
