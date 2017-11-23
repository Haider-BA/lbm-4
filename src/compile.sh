cflags="-c -Wall -march=native -std=f2008 -O3"
lflags="-Wall -march=native -std=f2008 -O3"

gfortran $cflags *.f95
gfortran $lflags *.o -o main