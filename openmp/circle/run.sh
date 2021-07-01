g++ -o circle -fopenmp -DSIZE=100000000 circle.cpp && timeout 60s time ./circle
