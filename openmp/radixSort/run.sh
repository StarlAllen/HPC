g++ -o sort_radix -fopenmp -DSIZE=1000000 sort_radix.cpp && timeout 60s time ./sort_radix
