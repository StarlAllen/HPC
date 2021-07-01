g++ -o scan_link -fopenmp -DSIZE=1000000 scan_link.cpp && timeout 60s time ./scan_link
