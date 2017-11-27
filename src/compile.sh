cflags="-c -Wall -march=native -std=f2008 -O3"
lflags="-Wall -march=native -std=f2008 -O3"

gfortran $cflags precision_mod.f95
gfortran $cflags math_mod.f95
gfortran $cflags utils_mod.f95
gfortran $cflags io_mod.f95
gfortran $cflags rootfinder_mod.f95
gfortran $cflags lagrange_mod.f95
gfortran $cflags boundary_condition_mod.f95
gfortran $cflags cell_mod.f95
gfortran $cflags sman_cell_mod.f95
gfortran $cflags abstract_lattice_mod.f95
gfortran $cflags lattice_mod.f95
gfortran $cflags periodic_lattice_mod.f95
gfortran $cflags sman_lattice_mod.f95
gfortran $cflags main.f95

gfortran $lflags *.o -o main

rm *.o
rm *.mod