
#define CUTOFF 1024

#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

// add your codes begin
void print(vector<int> vec){
    for(int i=0;i<vec.size();i++)
      printf("%d%c",vec[i],i< vec.size()-1? ' ':'\n');
}
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
//	for(i=0;i<SIZE;i++)
//		cout<<data[i]<<" ";
	//cout<<endl<<"num_threads:"<<num_threads<<endl;
}
void scan(vector<int> flag,vector<int>& data,vector<int> ptr, vector<int> tmpdata, vector<int> tmpptr) {
  bool done = false;
  while (!done ) {
    done = true;
    #pragma omp parallel for reduction(and:done)
    for (int i = 0; i < SIZE; i++) {
      if (ptr[i] >= 0) {
        tmpdata[i] = data[i] + data[ptr[i]];
        tmpptr[i] = ptr[ptr[i]];
        done = false;
      }
    }
    #pragma omp parallel for
    for (int i = 0; i < SIZE; i++) {
      data[i] = tmpdata[i];
      ptr[i] = tmpptr[i];
    }
  }
}

vector<int>  compact(vector<int> &data,vector<int>& psum,vector<int> &flag){
    vector<int> comp1(psum[SIZE-1],0);
    #pragma omp parallel for
    for(int i=0;i < SIZE;i++){
        if( flag[i] == 1 ){
            comp1[psum[i]-1] = data[i];
        }
    }
    return comp1;
}

void radix_sort(vector<int>& data){
  int bits = ceil(log2(SIZE*10));
  vector<int> flag(SIZE,0);
  vector<int> prev(SIZE,-1);
  vector<int> psum(SIZE,1);

  for(int i=0;i<bits;i++){
    //int count1=0;
    #pragma omp parallel for
    for(int j=0;j<SIZE;j++){ 
      psum[j] = 1-((data[j]>>i)&1);
      flag[j] = psum[j];
      prev[j] = j-1;
    }    
    //printf("%d flag:  ",i);print(flag);   
    //from exp of scan_vec
    inclusive_scan(psum);
    //scan(flag,psum,prev,psum,prev);
    int num1 = psum[SIZE-1];
    vector<int> tmpdata(SIZE); 
    //tmpdata = compact(data,psum,flag);
    #pragma omp parallel for
    for(int j=0;j< SIZE;j++){
        if(flag[j]==1){
          tmpdata[psum[j]-1] = data[j];
        }
        else
          tmpdata[num1 + j-psum[j]] = data[j];      
    }
    for(int i=0;i<SIZE;i++){
    	data[i] = tmpdata[i];
	} 
    //data = tmpdata;    
    //vector<int> comp0 = compact(data,psum,flag0);
    /*
    #pragma omp parallel for
    for(int j=0;j < SIZE;j++){
      flag[j] = 1-flag[j];
      psum[j] = flag[j];
      prev[j] = j-1;
    }
    */
    //scan(psum,prev,psum,prev);
    //vector<int> comp1 = compact(data,psum,flag0);
    //printf("%d data ",i);print(data);    
  }
}
// add your codes end
int main() {
  vector<int> data(SIZE);

  srand(SIZE);
  for (int i = 0; i < SIZE; i++) data[i] = rand() % (SIZE * 10);

  double t = omp_get_wtime();
  // add your codes begin
  radix_sort(data);
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  for (int i = 0; i < SIZE-1; i++) assert(data[i] <= data[i+1]);
}

