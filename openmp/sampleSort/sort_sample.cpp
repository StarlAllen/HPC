
#define CUTOFF 1024

#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

// add your codes begin
#include <cmath>
//compare function used for qsort
int cmp(const void *a,const void *b){
  return (*(int*)a-*(int*)b);
}
void sort_sample(vector<int> &data){  
  vector<int> splitters,selected_splitters;
  vector< vector< vector<int> >> bkt;
  int tid,i,j,k,idx,bkt_size,bkt_len,bk_idx,bkt_gap;

// split bin
#pragma omp parallel private(tid,i) shared(k)
{
    #pragma omp single
    {  
    k = omp_get_num_threads();
    //printf("num_threads:%d\n",k);
    bkt_size = ceil(SIZE/k);
    bkt_gap = (SIZE/(k*k));

    //Creating space for each thread to avoid dataracing(write) problem
    bkt.resize(k);
    for(i=0;i<k;i++){
     bkt[i].resize(k);
    } 
    // get splitters
    for(i=0;i<k*k;i++){
        if(i%k==0){
        bkt_len = (i+bkt_gap*k > SIZE ? SIZE-i:bkt_gap*k);
        qsort(&data[0]+i,bkt_len,sizeof(int),cmp);
        }
    }
    idx=0;
    for(int i=0;i<k*k;i++){
        if(i%k){
            splitters.push_back(data[i*bkt_gap]);
            idx++;
        }
    }    
    qsort(&splitters[0],splitters.size(),sizeof(int),cmp);   
    // choose splitters
    for(int i = 1;i < k; i++){
        selected_splitters.push_back(splitters[i*k-1]);
    }
    }//end of splitters
}
    //Data barrel
    #pragma omp parallel for private(tid,bk_idx)
    for(i = 0;i < SIZE;i++){
    bk_idx=0;
	tid = omp_get_thread_num();
        while(data[i] > selected_splitters[bk_idx] ){ 
            bk_idx++;
            if(bk_idx == k-1){
                break;
            }
        }
        bkt[tid][bk_idx].push_back(data[i]);//
    }
    //merge buckets
    for(i=1;i<k;i++)// threads
    	for(int j=0;j<k;j++) //bkt_idx
    	bkt[0][j].insert(bkt[0][j].end(),bkt[i][j].begin(),bkt[i][j].end());
    
    //qsort for each bucket
    #pragma omp parallel private(tid)
    {   
        tid = omp_get_thread_num();
	    qsort(&bkt[0][tid][0],bkt[0][tid].size(),sizeof(int),cmp);
    }
  data.clear();
  for(int i=0;i<k;i++){
      data.insert(data.end(),bkt[0][i].begin(),bkt[0][i].end());
  }
}
// add your codes end

int main() {
  vector<int> data(SIZE);

  srand(SIZE);
  for (int i = 0; i < SIZE; i++) data[i] = rand() % (SIZE * 10);

  double t = omp_get_wtime();
  // add your codes begin
  sort_sample(data);
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  for (int i = 0; i < SIZE-1; i++) assert(data[i] <= data[i+1]);
}

