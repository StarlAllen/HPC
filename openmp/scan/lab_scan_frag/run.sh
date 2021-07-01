g++ -o scan_frag -fopenmp -DSIZE=1000000 scan_frag.cpp && timeout 60s time ./scan_frag
