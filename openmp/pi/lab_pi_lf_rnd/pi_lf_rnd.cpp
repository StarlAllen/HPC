
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <vector>
using namespace std;

static unsigned long long MULTIPLIER  = 764261123;
static unsigned long long PMOD        = 2147483647;
static unsigned long long mult_n;

#define MAX_THREADS 128
static unsigned long long pseed[MAX_THREADS][4];//[4] to padding to cache line
                                                //size to avoid false sharing
// add your codes begin
long random_last = 0;
#pragma omp threadprivate(random_last)
double my_random()
{
       
	long random_next = (mult_n*random_last)%PMOD;
	random_last = random_next;
	return ((double)random_next/(double)PMOD);
}
double PI()
{
	long i,Ncirc=0;
	double pi,r=1.0;
	double x,y;
	int id;
#pragma omp parallel
{

#pragma omp single
{
	int nthreads =omp_get_num_threads();
	unsigned long long iseed = PMOD/MULTIPLIER;
	pseed[0][0] = iseed;
	mult_n = MULTIPLIER;
	for(int k=1;k<nthreads;++k)
	{
	iseed = (unsigned long long)((MULTIPLIER*iseed)%PMOD);
	pseed[k][0] = iseed;
	mult_n = (mult_n*MULTIPLIER)%PMOD;
	}
}       
    id = omp_get_thread_num();
	random_last = (unsigned long long) pseed[id][0];               

#pragma omp  for private(x,y) reduction(+:Ncirc)
	for(i=0;i<SIZE;i++)
	{
	x = my_random();y = my_random();	
    	if( x*x+y*y <= r*r)
		Ncirc++;
	}
}//end of  parallel
	pi = 4.0*((double)Ncirc/(double)SIZE);
	return pi;

}

// add your codes end
int main() {
  double pi;

  double t = omp_get_wtime();
  srand(20200218);
  // add your codes begin
  pi = PI();
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  printf("pi %.12f %.12f\n", pi, pi-M_PI);
}
