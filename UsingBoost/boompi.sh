mpicxx $1.cpp -o $1  -lmpi -lmpi_cxx -I/usr/include/mpich
mpirun -np $2 ./$1

