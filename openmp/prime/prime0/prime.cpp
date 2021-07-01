
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;


// add your codes begin
typedef  unsigned int unit32;
typedef unsigned long long int unit64;
#define maxn_num_threads 64
void par_prime(vector<long> prime){
  
  vector<int> p(SIZE+1,1);
  vector< vector<int>> p_n(maxn_num_threads);
  int num_threads,tid;
  #pragma omp parallel for
  for(unit64 i=2;i<=(unit64)sqrt(SIZE);i++){
    if(p[i]){
      for(unit64 j=i;i*j<=SIZE;j++)
        p[i*j] = 0;
    }
  }
  #pragma omp parallel
  {
    num_threads = omp_get_num_threads();
    tid = omp_get_thread_num();
    //ceil(SIZE/num_threads)
    unit64 st = tid*(SIZE/num_threads);
    unit64 ed = (tid==num_threads-1 ? SIZE-1:(tid+1)*(SIZE/num_threads));
    for(unit64 i=st;i<ed;i++){
      if(p[i]) p_n.push_back(i);
    }
  }
  for(int k=0;k<num_threads;k++){
    prime.insert(prime.end(),p_n[k].start(),p_n[k].end());
  }
  
}
// add your codes end


int main() {
  vector<long> prime;

  double t = omp_get_wtime();
  // add your codes begin
  par_prime(prime)
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %ld\n", t, long(SIZE));

  printf("prime");
  sort(prime.begin(), prime.end());
  printf("\nsize %ld\n", prime.size());
}

