#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <iostream>
using namespace std;


// add your codes begin
#include <cmath>
void print(vector<int> vec){
     for(int i=0;i<vec.size();i++){
     	printf("%d ",vec[i]);
     }
     printf("\n");
     return;
}

void scan_link(vector<int> &data,vector<int> data_old,vector<int> &prev,vector<int> pre_old){

	//for(int j=0;j<ceil(log2(SIZE));j++){
	bool done = false;
	while(!done){
		done = true;
		#pragma omp parallel
		{
		#pragma omp for
     		for(int i=0;i<SIZE;i++){
     			if(prev[i]!=-1){
     				done = false;
					data_old[i] += data[prev[i]];
					pre_old[i] = prev[prev[i]];
        		}
     		}
     	#pragma omp for
			for(int i = 0;i < SIZE;i++){			
     			prev[i] = prev_old[i];
				data[i] = data_old[i];     
			}
   		}
    }
}
// add your codes end


int main() {
  vector<int> data(SIZE, -1);
  vector<int> prev(SIZE, -1);
  vector<int> next(SIZE, -1);
  vector<int> test(SIZE, -1);

  srand(SIZE);
  { int tmp = -1;
    for (int i = 0; i < SIZE/10; i++) {
      int idx = rand() % SIZE;
      while (data[idx] >= 0) idx = (idx + 1) % SIZE;
      if (i > 0) {
        data[idx] = 1;
        prev[idx] = tmp;
        next[tmp] = idx;
      } else {
        data[idx] = 0;
      }
      test[idx] = i;
      tmp = idx;
    }
  }

  double t = omp_get_wtime();
  // add your codes begin
  scan_link(data,data,prev,pre);
//  printf("data next prev test\n");
//  print(data);
//  print(next);
//  print(prev);
//  print(test);
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  for (int i = 0; i < SIZE; i++) assert(data[i] == test[i]);
}

