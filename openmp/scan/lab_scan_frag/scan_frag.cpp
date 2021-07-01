#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

// add your codes begin
template<class T>
void print(vector<T> vec){
    for(int i=0;i<vec.size();i++){
      cout<<vec[i]<<" ";
            // printf("%d%c",data[i],i==data.size()-1? '\n':' '));
    }
    cout<<endl;
    //printf("\n");
}
void scan_frag(vector<int> &data,vector<bool> &flag,vector<int> &pre , vector<bool> flag_init )
{
        int d,k,alen=data.size();
        data.resize(pow(2,log2(alen-1)+1),0);
        flag.resize(pow(2,log2(alen-1)+1),1);
        flag_init.resize(pow(2,log2(alen-1)+1),1);
                //printf("flag_init:");
        int n = data.size();
        //up-sweep
        //printf("\n%d\n",(int)log2(n)-1);
        for(d=0;d<=((int)log2(n)-1);d++){
               #pragma omp parallel
                {
                #pragma omp for
                for(k=0;k<= n-(int)pow(2,d+1);k+=(int)pow(2,d+1) ){
                //      if(k+pow(2,d+1)-1>=n) continue;
                        if(flag[k+pow(2,d+1)-1]==0)
                             data[k+pow(2,d+1)-1] = data[k+pow(2,d)-1]+data[k+pow(2,d+1)-1];
                        flag[k+pow(2,d+1)-1] = flag[k+pow(2,d)-1]||flag[k+pow(2,d+1)-1];

                }
                }
//           printf("%d:\n",d);
//           print(flag);
//           print(data);
        }
//      printf("up-sweep end\ndown-sweep start\n");
        //print(data); 
        //down-sweep
        data[n-1] = 0;
        for(d=((int)log2(n)-1);d >= 0;d--){
           #pragma omp parallel
                {
                #pragma omp for
                for(k=0;k <=(int)( n-pow(2,d+1));k+=(int)pow(2,d+1)){
                     int tmp=data[k+pow(2,d)-1];
                     data[k+pow(2,d)-1] = data[k+pow(2,d+1)-1];
                     if( flag_init[k+pow(2,d)] == 1)
                            data[k+pow(2,d+1)-1] = 0;
                     else if( flag[k+pow(2,d)-1] == 1)
                            data[k+pow(2,d+1)-1] = tmp;
                     else
                            data[k+pow(2,d+1)-1] = tmp + data[k+pow(2,d+1)-1];
                     flag[k+pow(2,d)-1] = 0;
                }
              }
//          printf("%d:\t",d);
//          print(data);
        }
     //  print(data);
       #pragma omp parallel for
           for(int i = 1;i< alen;i++){
                   if(!flag_init[i])
                      pre[i] += data[i];
           }
}
// add your codes end

int main() {
  vector<int> data(SIZE, 1);
  vector<bool> flag(SIZE, false);
  vector<int> test(SIZE);

  srand(SIZE);
  data[0] = 0; flag[0] = true;
  for (int i = 0; i < flag.size() / 12; i++) {
    int index = rand() % flag.size();
    data[index] = 0; flag[index] = true;
  }
  
  for (int i = 0; i < data.size(); i++) test[i] = flag[i] ? data[i] : test[i-1] + data[i];

  double t = omp_get_wtime();
  // add your codes begin
  
  vector<int> pre(data);
  vector<bool> flag_init(flag);
//  printf("data,flag,test,data\n");
//  print(data);
//  print(flag);
//  print(test);
   scan_frag(data,flag,pre,flag_init);
//  print(data);
//   print(test);
//   print(pre);
  // add your codes end
  t = omp_get_wtime() - t;
  printf("time %f %d\n", t, SIZE);

  for (int i = 0; i < SIZE; i++) assert(pre[i] == test[i]);
       // if(pre[i]!=test[i])  printf("%d\n",i);
}  

