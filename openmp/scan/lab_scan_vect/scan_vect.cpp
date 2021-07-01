#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <iostream>
using namespace std;

// add your codes begin
/*
void exclusive_scan_iterative(int* data, int size)
{
        // upsweep phase.
        for (int twod = 1; twod < size; twod*=2)
        {
                int twod1 = twod*2;
                //parallel_for (int i = 0; i < size; i += twod1)
                #pragma omp parallel for{
                        for(int i=0;i<size;i+= twod1)
                                data[i+twod1-1] += data[i+twod-1];
                }
        }
        data[N-1] = 0;
        // downsweep phase.
        for (int twod = size/2; twod >= 1; twod /= 2)
        {
                        // downsweep phase.
        for (int twod = size/2; twod >= 1; twod /= 2)
        {
                int twod1 = twod*2;
                //parallel_for (int i = 0; i < size; i += twod1)
                #pragma omp parallel for
                for(int i=0;i<size; i+= twod1)
                {
                        int t = data[i+twod-1];
                        data[i+twod-1] = data[i+twod1-1];
                        // change twod1 below to twod to reverse prefix sum.
                        data[i+twod1-1] += t;
                }
        }
}
*/
void inclusive_scan(vector<int> &data){

        int tid,i,j,seq_num=0,num_threads=0;
        //step  1
        #pragma omp parallel private(tid,j) shared (data)
        {
                tid = omp_get_thread_num();
                num_threads = omp_get_num_threads();
                seq_num = SIZE/num_threads;

                for(j = 0;j<seq_num -1;j++){
                        data[tid*seq_num+j+1] += data[tid*seq_num+j];
                   }
        }

        // step 2
        for(i=2;i<num_threads+1;i++){
                data[i*seq_num-1] += data[(i-1)*seq_num-1];
        }
        //step 3

        #pragma omp parallel private(tid,j) shared(data)
        {
                tid = omp_get_thread_num();
                if(tid!=0){
                        for(j=0;j<seq_num-1;j++){
                                data[tid*seq_num+j] += data[tid*seq_num-1];
                        }
                }
        }
        //step 4
       if(SIZE%num_threads){
                for(i=num_threads*seq_num;i<SIZE;i++)
                    data[i] += data[i-1];
        }
//      for(i=0;i<SIZE;i++)
//              cout<<data[i]<<" ";
        //cout<<endl<<"num_threads:"<<num_threads<<endl;
}
// add your codes end

int main() {
  vector<int> data(SIZE, 1);
  data[0] = 0;

  double t = omp_get_wtime();
  // add your codes begin
  inclusive_scan(data);
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  for (int i = 0; i < SIZE; i++) assert(data[i] == i);
}
