g++ -o prime -fopenmp -DSIZE=1000000 prime.cpp && timeout 60s time ./prime
