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
/*
void scan_link(vector<int> &data,vector<int> &prev,vector<int> &next){

	vector<int> data_old(data);
  vector<int> pre_old(prev);
  
  for(int j=0;j<ceil(log2(SIZE));j++){
	#pragma omp parallel
	{
		#pragma omp for
     		for(int i=0;i<SIZE;i++){
     			if(prev[i]!=-1){
			data_old[i] += data[prev[i]];
			pre_old[i] = prev[prev[i]];
        		}
     		}
    prev.assign(pre_old.begin(),pre_old.end());     
		data.assign(data_old.begin(),data_old.end());
   	}
  }
}
*/

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

}      

void scan_circle_link(vector<int> &data,vector<int> &prev,vector<int> data_old,vector<int> prev_old,vector<int> &max_idx,vector<int> max_idx1,vector<int> flag){//,vector<int> &next){
	
 /*
  double time = omp_get_wtime();
 
  vector<int> data_old(SIZE);
  //vector<int> tmpprev(prev);
  vector<int> prev_old(SIZE);
  vector<int> flag(SIZE);
  //vector<int> flag_old(SIZE,0);
  //vector<int> flag1(SIZE,0);
  vector<int> max_idx(SIZE);// 作为分组依据
  vector<int> max_idx1(SIZE);
  #pragma omp parallel for
  for(int i=0;i<SIZE;i++){
    data_old[i] = data[i];
    prev_old[i] = prev[i];
    if(data[i] != -1 )//&& i < next[i])
    {
      max_idx[i] = i; 
      max_idx1[i] = i; 
    }else{
      max_idx[i] = -1;
      max_idx1[i] = -1;
    }
  }
  time = omp_get_wtime()-time;
  printf("%lf\n",time);
  */

  //nlogn  maxscan
//double time = omp_get_wtime();
  //max-scan on index
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
//printf("max_scan:%lf\n",omp_get_wtime()-time);  
//time = omp_get_wtime(); 
  int cnt = 0;
  //#pragma omp parallel for reduction(+:cnt)
  //这里我们要根据flag得到compact下表，表比较稀疏而且会出现数据写冲突，因此采用串行遍历即可 且只占程序运行时间的一小部分
  for(int i=0;i<SIZE;i++){
    if(data[i]!= -1 && max_idx[i] == i){
      //tmpprev[next[i]] = -1; //转化成单链表的scan问题
      flag[i] = cnt;
      cnt++;
      /*
      int j = i;
      int nums = 1;
      while(next[j] != i){
        nums++;
        j = next[j];
      }
      data_old[i] = nums;
      */
    }
  }
//printf("cnt:%lf\n",omp_get_wtime()-time);
  //scan(data_old,tmpprev);
//printf("addscan:%lf\n",omp_get_wtime()-time);
//time = omp_get_wtime();
  int res[max_threads_nums][cnt];
  memset(res,0,sizeof(res));
  #pragma omp parallel for
  for(int i=0;i<SIZE;i++){
    int tid = omp_get_thread_num();
    if(data_old[i] != -1 ){
      //printf("%d ",data[i]);
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
//printf("write to data:%lf\n",omp_get_wtime()-time);
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

