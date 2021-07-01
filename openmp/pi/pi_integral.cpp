#include <omp.h>
#include <math.h>
#include <stdio.h>
using namespace std;

// add your codes begin
#include<cstring>
#define MAX_THREADS 64
#define num_threads 64
#define PAD 4

double sum[MAX_THREADS][PAD];//[4];
//double sum[MAX_THREADS];
double step,pi,x;
double PI()
{
    int i,id,nthreads;
    step = 1.0/(double)SIZE;
    omp_set_num_threads(num_threads);
    #pragma omp parallel private(i,id,x)
    {
        id = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        if(id==0) nthreads = nthrds;
        
	for(i=id,sum[id][0]=0.0;i<SIZE;i+= nthrds){
            x = (i+0.5)*step;
            sum[id][0] += 4.0/(1.0+x*x);
        }/*
	for(i=id,sum[id]=0.0;i<SIZE;i+= nthrds){
            x = (i+0.5)*step;
            sum[id] += 4.0/(1.0+x*x);
        }*/
    }
    pi=0.0;
    for(i=0;i<nthreads;i++)
        pi += sum[i][0];
       //pi += sum[i];
    pi*= step;
    return pi;
}

// add your codes end


int main() {
  double pi;

  double t = omp_get_wtime();
  // add your codes begin
  pi = PI();
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  printf("pi %.12f %.12f\n", pi, pi-M_PI);
}

