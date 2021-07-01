g++ -o sort_sample -fopenmp -DSIZE=100000000 sort_sample.cpp && timeout 60s time ./sort_sample
