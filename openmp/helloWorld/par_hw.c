
#include "stdio.h"
#include "omp.h"

int main(){
   int cnt=0;
   #pragma omp parallel
   {
       int ID = omp_get_thread_num();
	printf("hello((%d)",ID);
        printf("world(%d)\n",ID);
	cnt = omp_get_num_threads();
    }//return 0;
   printf("serial_NUM_threads:%d\n",omp_get_num_threads());
   printf("parallel_num_threads:%d\n",cnt);
}

/*
int main(){
    int count = 0;
    //int ID =omp_get_thread_num();
    #pragma omp parallel  num_threads(2)
    {
        for(int i=0;i<10000;i++){
            count++;
        }
    }
    printf("%d(%d)\n",count,omp_get_thread_num());
}
*/
