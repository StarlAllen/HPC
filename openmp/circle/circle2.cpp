
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <algorithm>
using namespace std;


// add your codes begin
#include <cmath>
#include <cstring>
#define max_threads_nums 128

//scan  linked list in  array pool 
void scan(vector<int>& data, vector<int> ptr) {
  vector<int> tmpdata(data);
  vector<int> tmpptr(ptr);
  bool done = false;
  while (!done) {
    done = true;
    #pragma omp parallel for reduction(and:done)
    for (int i = 0; i < SIZE; i++) {
      if (ptr[i] >=0) {
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

void scan_circle_link(vector<int> &data,vector<int> &prev,vector<int> data_old,vector<int> prev_old,vector<int> &max_idx,vector<int> max_idx1,vector<int> flag){//,vector<int> &next){

/*
I follow three steps to calculate circle nums and their sizes:
1.max_scan  find the max_index of each list(circle) and sign them.
2.map  count the max_index nums as cnt and set mapping for different elements from  max_index  to array(0...cnt)
3.parallel  add&&write to array
*/

//1.max-scan on index
  bool flg0 = false;
  while(!flg0){
    flg0 = true; 
    #pragma omp parallel for reduction(and:flg0)
    for(int j=0;j<SIZE;j++){
      if(data[j] != -1){
        if(max_idx[j] < max_idx[prev[j]]){
            flg0 = false;
            max_idx1[j] = max_idx[prev[j]];
            prev_old[j] = prev[prev[j]];
        }
      }
    }
    #pragma omp parallel for
    for(int i=0;i<SIZE;i++){
        max_idx[i] = max_idx1[i];
        prev[i] = prev_old[i];
    }
  }
//2.mapping from max_idx to array_idx
  int cnt = 0;
  //这里我们要根据flag得到compact下表，表比较稀疏而且会出现数据写冲突，因此采用串行遍历即可 且只占程序运行时间的一小部分
  for(int i=0;i<SIZE;i++){
    if(data[i]!= -1 && max_idx[i] == i){
      flag[i] = cnt;
      cnt++;
    }
  }
//3.add && write to array 
  int res[max_threads_nums][cnt];
  memset(res,0,sizeof(res));
  #pragma omp parallel for
  for(int i=0;i<SIZE;i++){
    int tid = omp_get_thread_num();
    if(data_old[i] != -1 ){
      res[tid][flag[max_idx[i]]] += data_old[i];
    }
  }
  data.resize(cnt);
  for(int j=0;j<cnt;j++){
    data[j] = 0;
    for(int i=0;i<max_threads_nums;i++){
      data[j] += res[i][j];
    } 
  } 
}
// add your codes end

int main() {
  vector<int> data(SIZE, -1);
  vector<int> prev(SIZE, -1);
  vector<int> next(SIZE, -1);
  vector<int> test;

  srand(20200218);
  { int empty = SIZE / 10;
    int head, tail, tmp;
    while (empty > 0) {
      int size = rand() % empty + 1;
      for (int i = 0; i < size; i++) {
        int idx = rand() % SIZE;
        while (data[idx] >= 0) idx = (idx + 1) % SIZE;

        data[idx] = 1;
        if (i == 0) {
          head = idx;
          tail = idx;
        } else if (i == size-1) {
          prev[idx] = tmp; next[tmp] = idx;
          tail = idx;
        } else {
          prev[idx] = tmp; next[tmp] = idx;
        }

        tmp = idx;
      }
      prev[head] = tail; next[tail] = head;
      test.push_back(size);
      empty -= size;
    }
    sort(test.begin(), test.end());
  }

  double t = omp_get_wtime();
  // add your codes begin
//scan_circle_link(vector<int> &data,vector<int> &prev,vector<int> data_old,vector<int> prev_old,vector<int> &max_idx,vector<int> max_idx1,vector<int> flag){//,vector<int> &next){
//double time = omp_get_wtime();
  vector<int> max_idx(SIZE);
  vector<int> flag(SIZE);
  #pragma omp parallel for
  for(int i=0;i<SIZE;i++)
     if(data[i]!=-1){
       max_idx[i] = i;
     }else{
       max_idx[i] = -1;
     }
//printf("%lf\n",omp_get_wtime()-time);
  scan_circle_link(data,prev,data,prev,max_idx,max_idx,flag);
  // printf("\n test:");
  // for (int i = 0; i < test.size(); i++) printf(" %d", test[i]);
  //   printf("\nsize %d\n", int(test.size()));  
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  sort(data.begin(), data.end());
  printf("circle");
  for (int i = 0; i < data.size(); i++) printf(" %d", data[i]);
  printf("\nsize %d\n", int(data.size()));
  assert(data.size() == test.size());
  for (int i = 0; i < test.size(); i++) assert(data[i] == test[i]);
}

