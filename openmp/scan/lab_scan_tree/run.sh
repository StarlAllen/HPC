g++ -o scan_tree -fopenmp -DSIZE=1000000 scan_tree.cpp && timeout 60s time ./scan_tree
