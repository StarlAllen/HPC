#include <omp.h>
#include <math.h>
#include <stdio.h>
using namespace std;

// add your codes begin
#include<cstring>
static unsigned long long MULTIPLIER  = 764261123;
static unsigned long long PMOD        = 2147483647;
static unsigned long long mult_n;

#define MAX_THREADS 128
static unsigned long long pseed[MAX_THREADS][4];
//static unsigned long long pseed[MAX_THREADS];                       
unsigned long long random_last = 0;
#pragma omp threadprivate(random_last)

double my_random()
{
       
	long random_next = (mult_n*random_last)%PMOD;
	random_last = random_next;
	return ((double)random_next/(double)PMOD);
}
double PI()
{
  long Ncirc = 0;
  long cir[MAX_THREADS][4];//[4] to padd to cache line
//  long cir[MAX_THREADS];
  double pi,r=1.0,x,y;
  int i,nthreads;
  memset(cir,0,sizeof(cir));
#pragma omp parallel
{

#pragma omp single
{
	nthreads =omp_get_num_threads();
	unsigned long long iseed = PMOD/MULTIPLIER;
//	pseed[0] = iseed;	
	pseed[0][0] = iseed;
	mult_n = MULTIPLIER;
	for(int k=1;k<nthreads;++k)
	{
	iseed = (unsigned long long)((MULTIPLIER*iseed)%PMOD);
//	pseed[k] = iseed;
	pseed[k][0] = iseed;
	mult_n = (mult_n*MULTIPLIER)%PMOD;
	}
}
int id =  omp_get_thread_num();
random_last = (unsigned long long) pseed[id];

#pragma omp  for private(x,y)
	for(i=0;i < SIZE ;i++)
	{    
    x = my_random();y = my_random();	
    if( x*x+y*y <= r*r)
		{
      cir[id][0]+=1;
//`        cir[id] += 1;
		}
	}
}
//end of  parallel
  for(i=0;i<nthreads;i++)
    {
   Ncirc += cir[i][0];
//     Ncirc += cir[i];
    //printf("%d:%ld\n",i,cir[i][0]);
    }
	pi = 4.0*((double)Ncirc/(double)SIZE);
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

