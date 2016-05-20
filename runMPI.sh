mpicxx $1.cpp -o ./build/$1
mpirun -np $2 ./build/$1
