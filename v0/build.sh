# compile objects
gfortran -c -O3 -march=native -Wall precision_mod.f95 utils_mod.f95 rootfinder_mod.f95 newunit_mod.f95 io_mod.f95 lbm_mod.f95 lbm1d_mod.f95 

# compile driver
gfortran -O3 -Wall -march=native *.o heat_equation.f95 -o heat_equation