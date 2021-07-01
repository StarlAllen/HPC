#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <vector>
using namespace std;

static long MULTIPLIER  = 1366;
static long ADDEND      = 150889;
static long PMOD        = 714025;

// add your codes begin
long random_last=0;
#pragma omp threadprivate(random_last) //for threadsafe
double MY_random()
{
	long random_next;
	random_next = (MULTIPLIER *random_last+ADDEND)%PMOD;
	random_last = random_next;
	return ((double)random_next/(double)PMOD);
}
double PI()
{
	long i,Ncirc=0;
	double x,y,pi,r=1.0;
	//seed(0,-r,r);
#pragma omp parallel for private(x,y) reduction(+:Ncirc)
	for(i=0;i<SIZE;i++)	
	{
		//x=random();y=random();
		x = MY_random();y = MY_random();
		if(x*x + y*y <=r*r)
			Ncirc++;
	}
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

