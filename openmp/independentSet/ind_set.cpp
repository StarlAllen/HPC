#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;


// add your codes begin
// add your codes end

static long MULTIPLIER  = 1366;
static long ADDEND      = 150889;
static long PMOD        = 714025;

long rnd(long& seed) {
  seed = (MULTIPLIER  * seed + ADDEND) % PMOD;
  return seed;
}

void scan(vector<int>& data, vector<int> ptr, vector<int> tmpdata, vector<int> tmpptr) {
  bool done = false;
  while (!done) {
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


int main() {
  vector<int> data(SIZE, -1);
  vector<int> prev(SIZE, -1);
  vector<int> next(SIZE, -1);
  vector<int> test(SIZE, -1);

  vector<int> chkidx(SIZE, -1);
  vector<int> chkdata(SIZE, -1);
  vector<int> chkprev(SIZE, -1);
  vector<int> chknext(SIZE, -1);

  srand(20200218);
  int size = SIZE / 10;
  long seed[omp_get_max_threads()][8];
  for (int i = 0; i < omp_get_max_threads(); i++) seed[i][0] = rand();
  { int tmp = -1;
    for (int i = 0; i < size; i++) {
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
  int idx = 0, count = size;
  for (; count > int(size / log(size)); idx++) {
    printf("iter %d\n", idx);
    vector<int> flag(SIZE, 0);
    #pragma omp parallel for
    for (int i = 0; i < SIZE; i++) {
      if ((data[i] >= 0) && (prev[i] >= 0) && (next[i] >= 0)) {
        flag[i] = rnd(seed[omp_get_thread_num()][0]) % 2;
      }
    }
    #pragma omp parallel for
    for (int i = 0; i < SIZE; i++) {
      if ((flag[i] != 0) && (next[i] >= 0) && (flag[next[i]] != 0)) {
        flag[i] = 0;
      }
    }

    // add your codes begin
    #pragma omp parallel for reduction(-:count)
	for(int i=0;i<SIZE;i++){
    	if(flag[i]==1){
    		if(prev[i]>=0 ) next[prev[i]] = next[i];
			if(next[i]>=0){
				data[next[i]] += data[i];
				prev[next[i]] =  prev[i];
			}
			chkdata[i] = data[i];
			chkprev[i] = prev[i];
			chknext[i] = next[i];
			data[i] = -1;
			chkidx[i] = idx;
			count--;  
		}
	} 
    // add your codes end
  }

  int totaldata = 0, countdata = 0;
  #pragma omp parallel for reduction(+:totaldata) reduction(+:countdata)
  for (int i = 0; i < SIZE; i++) {
    if (data[i] >= 0) {
      totaldata += data[i];
      countdata += 1;
    }
  }
  assert(totaldata == size-1);
  assert(countdata < size);
  printf("delete %.1f%%\n", 100.0 * (size - countdata) / size);

  scan(data, prev, data, prev);

  idx = idx - 1;
  for (; idx >= 0; idx--) {
    printf("iter %d\n", idx);
    // add your codes begin
     
	#pragma omp parallel for
    for(int i=0;i<SIZE;i++){
    	if(chkidx[i]== idx){
    		next[chkprev[i]] = i;
    		prev[chknext[i]] = i;
    		data[i] = data[next[i]]-chkdata[i]; 
		}
	}
	// add your codes end
  }
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  for (int i = 0; i < SIZE; i++) assert(data[i] == test[i]);
}

