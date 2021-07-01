g++ -o pi_lcg_rnd -fopenmp -DSIZE=1000000 pi_lcg_rnd.cpp && timeout 60s time ./pi_lcg_rnd
