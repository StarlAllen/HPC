
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

// add your codes begin
#define maxn_num_threads 64
vector<int> p(SIZE+1,1);
vector< vector<long> > p_n(maxn_num_threads);
int num_threads,tid;

void par_prime(vector<long>& prime){
//埃氏筛法
p[0]=p[1]=0;
//double t = omp_get_wtime();
/*
// this parallel method is slower than another
#pragma omp parallel
  {
  #pragma omp  for
  for(long i=2;i<=long(sqrt(SIZE));i++){
    if(p[i]){	    
      for(long j=i*i;j<=SIZE;j+=i)
       { 
        p[j] = 0;
       }
      }
    }
  }
*/
  for(long i=2;i<=long(sqrt(SIZE));i++){
    if(p[i]){
	 #pragma omp  parallel for   
      for(long j=i*i;j<=SIZE;j+=i)
       { 
        p[j] = 0;
       }
      }
    }
 //Euler筛法 don't know how to parallel it...
 /*
 for(long i=2;i<SIZE;i++)
 {
   if(p[i]) prime.push_back(i);
   for(long j=0;j<prime.size();j++){
     if(i*prime[j]>SIZE) break; 
        p[i*prime[j]] = 0;
     if(i%prime[j]==0)
        break;
   }
 }
 */
//printf("%f\n",omp_get_wtime()-t);
//t = omp_get_wtime();
if(SIZE < 1e5) 
  //serial
  for(long i=0;i<SIZE;i++){
    if(p[i]) prime.push_back(i);
  }
else{
#pragma omp parallel
 {
	#pragma omp single 
	 num_threads = omp_get_num_threads();
        #pragma omp  for private(tid)
  for(long i=0;i<SIZE;i++){
	if(p[i]){
  	tid = omp_get_thread_num();
	  //printf("tid:%d\n",tid);
	  p_n[tid].push_back(i);
	}
  }
 }
  for(int k=0;k<num_threads;k++){
	 prime.insert(prime.end(),p_n[k].begin(),p_n[k].end());
  }
}  
  //printf("%f\n",omp_get_wtime()-t);
}
  // add your codes end
int main() {
  vector<long> prime;

  double t = omp_get_wtime();
  // add your codes begin
  par_prime(prime);
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %ld\n", t, long(SIZE));

  printf("prime");
  sort(prime.begin(), prime.end());
  printf("\nsize %ld\n", prime.size());
}

