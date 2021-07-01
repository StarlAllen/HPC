g++ -o scan_vect -fopenmp -DSIZE=1000000 scan_vect.cpp && timeout 60s time ./scan_vect
