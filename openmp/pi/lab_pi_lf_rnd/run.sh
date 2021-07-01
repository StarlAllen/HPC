g++ -o pi_lf_rnd -fopenmp -DSIZE=1000000 pi_lf_rnd.cpp && timeout 60s time ./pi_lf_rnd
