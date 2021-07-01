#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

// add your codes begin
void print(vector<int> vec){
	for(int i=0;i<vec.size();i++){
	    printf("%d ",vec[i]);
	}
	printf("\n");	
}
void scan_tree(vector<int> &data,vector<int> data1,vector<int> &pare,vector<int> pare1){
	
	bool flag = false;
	while(!flag){
		{
			flag = true;
			#pragma omp parallel  for
     		for(int i=0;i<SIZE;i++){
     			if(pare[i]!=-1){
     			flag = false;
				data1[i] += data[pare[i]];
				pare1[i] = pare[pare[i]];
        		}
     		}
     		#pragma omp  parallel for
     		for(int i=0;i<SIZE;i++){
     			pare[i] = pare1[i];
     			data[i] = data1[i];
			 }
 		}
 	}
}
// add your codes end

int main() {
  vector<int> data(SIZE, -1);
  vector<int> pare(SIZE, -1);
  vector<int> test(SIZE, -1);

  srand(SIZE);
  { vector<int> tmp;
    for (int i = 0; i < SIZE/10; i++) {
      int idx = rand() % SIZE;
      while (data[idx] >= 0) idx = (idx + 1) % SIZE;
      if (i > 0) {
        data[idx] = 1;
        pare[idx] = tmp[rand() % tmp.size()];
        test[idx] = test[pare[idx]] + data[idx];
      } else {
        data[idx] = 0;
        test[idx] = data[idx];
      }
      tmp.push_back(idx);
    }
  }

  double t = omp_get_wtime();
  // add your codes begin
  scan_tree(data,data,pare,pare);
//  printf("data  pare  test\n");
//  print(data);
//  print(pare);
//  print(test);
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  for (int i = 0; i < SIZE; i++) assert(data[i] == test[i]);
}

