
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <iostream>

#include "mpi.h"


using namespace std;

void randMat(int rows, int cols, float *&Mat) {
  Mat = new float[rows * cols];
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Mat[i * cols + j] = 1.0;
}

void openmp_sgemm(int m, int n, int k, float *&leftMat, float *&rightMat,
                  float *&resultMat) {
  // rightMat is transposed
#pragma omp parallel for
  for (int row = 0; row < m; row++) {
    for (int col = 0; col < k; col++) {
      resultMat[row * k + col] = 0.0;
      for (int i = 0; i < n; i++) {
        resultMat[row * k + col] +=
            leftMat[row * n + i] * rightMat[col * n + i];
      }
    }
  }
  return;
}


void mpi_sgemm(int m, int n, int k, float *&leftMat, float *&rightMat,
               float *&resultMat, int rank, int worldsize) {

  int rowBlock = sqrt(worldsize);
  if (rowBlock * rowBlock > worldsize)
    rowBlock -= 1;
  int colBlock = rowBlock;

  int rowStride = m / rowBlock;
  int colStride = k / colBlock;

  worldsize = rowBlock * colBlock; // we abandom some processes.
  // so best set process to a square number.
  
  float *res;
//add code begin

  if (rank == 0) {
    float *buf = new float[k * n];
    // transpose right Mat
    for (int r = 0; r < n; r++) {
      for (int c = 0; c < k; c++) {
        buf[c * n + r] = rightMat[r * k + c];
      }
    }

    for(int r = 0; r < k; r++) {
      for (int c = 0; c < n; c++) {
        rightMat[r * n + c] = buf[r * n + c];
      }
    }
    delete buf;

    res = new float[rowStride*colStride];
    int sended_rows = 0,sended_cols = 0;
    for(int i=0;i<rowBlock;i++){
      for(int j=0;j<colBlock;j++){
        int pid = i*rowBlock + j;
        
	    rowStride = ((pid / colBlock) == rowBlock - 1)
                    ? m - (rowBlock - 1) * (m / rowBlock)
                    : m / rowBlock;
        colStride = ((pid % colBlock) == colBlock - 1)
                    ? k - (colBlock - 1) * (k / colBlock)
                    : k / colBlock;
        
        if(pid > 0){  //只有rank==0进程负责分发数据		        
	        MPI_Send(&leftMat[sended_rows*n+0],rowStride*n,MPI_FLOAT,pid,0,MPI_COMM_WORLD);
            MPI_Send(&rightMat[sended_cols*n+0],colStride*n,MPI_FLOAT,pid,1,MPI_COMM_WORLD);
        }
        sended_cols += colStride;
      }
      sended_rows += rowStride;
      sended_cols = 0;     
    }
  }
    
  if(rank > 0 && rank < worldsize){
	rowStride = ((rank / colBlock) == rowBlock - 1)
                    ? m - (rowBlock - 1) * (m / rowBlock)
                    : m / rowBlock;
    colStride = ((rank % colBlock) == colBlock - 1)
                    ? k - (colBlock - 1) * (k / colBlock)
                    : k / colBlock;    
    
    leftMat = new float[n*rowStride];
    rightMat = new float[n*colStride];
    res = new float[rowStride*colStride];
	    
    MPI_Recv(leftMat,n*rowStride,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(rightMat,n*colStride,MPI_FLOAT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
    }
  
//add code end
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank < worldsize) {
    rowStride = ((rank / colBlock) == rowBlock - 1)
                    ? m - (rowBlock - 1) * (m / rowBlock)
                    : m / rowBlock;
    colStride = ((rank % colBlock) == colBlock - 1)
                    ? k - (colBlock - 1) * (k / colBlock)
                    : k / colBlock;
    openmp_sgemm(rowStride, n, colStride, leftMat, rightMat, res);
	}

  MPI_Barrier(MPI_COMM_WORLD);

//add code begin
 //各个进程返回计算结果并返回给主进程
  if(rank > 0 && rank < worldsize){
        rowStride = ((rank / colBlock) == rowBlock - 1)
                    ? m - (rowBlock - 1) * (m / rowBlock)
                    : m / rowBlock;
        colStride = ((rank % colBlock) == colBlock - 1)
                    ? k - (colBlock - 1) * (k / colBlock)
                    : k / colBlock;

        // printf("rank:%d's res matrix\n",rank);
        // print(rowStride,colStride,res);
        MPI_Send(res,rowStride*colStride,MPI_FLOAT,0,1,MPI_COMM_WORLD);
        //MPI_Isend(res,rowStride*colStride,MPI_FLOAT,0,1,MPI_COMM_WORLD,&recvRequest[rank-1]);        
    }
  if(rank == 0){      
      int recv_row=0,recv_col=0;
      for(int i=0;i<rowBlock;i++){
        for(int j=0;j<colBlock;j++){
          int pid = i*rowBlock+j;          
          
          rowStride = ((pid / colBlock) == rowBlock - 1)
                    ? m - (rowBlock - 1) * (m / rowBlock)
                    : m / rowBlock;
          colStride = ((pid % colBlock) == colBlock - 1)
                    ? k - (colBlock - 1) * (k / colBlock)
                    : k / colBlock;          
          if(pid == 0){
              #pragma omp parallel for
			  for(int row=0;row<rowStride;row++){
                for(int col=0;col<colStride;col++){
                  resultMat[row*n+col] = res[row*colStride+col];
                }
              }
          }
          else{
            float *tmp = new float[colStride*rowStride];
            MPI_Recv(tmp,colStride*rowStride,MPI_FLOAT,pid,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            #pragma omp parallel for 
			for(int row=0;row<rowStride;row++){
              for(int col=0;col<colStride;col++){
                resultMat[n*(recv_row+row)+(recv_col+col)] = tmp[row*colStride+col];
              }
            }
            delete tmp;
          } 
          recv_col += colStride; 
        }
        recv_col = 0;
        recv_row += rowStride;          
      }
    }
//add code end
  MPI_Barrier(MPI_COMM_WORLD);
  return;
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    cout << "Usage: " << argv[0] << " M N K use-blas\n";
    exit(-1);
  }

  int rank;
  int worldSize;
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int k = atoi(argv[3]);
  //int blas = atoi(argv[4]);

  float *leftMat, *rightMat, *resMat;

  struct timeval start, stop;
  if (rank == 0) {
    randMat(m, n, leftMat);
    randMat(n, k, rightMat);
    randMat(m, k, resMat);
  }
  gettimeofday(&start, NULL);

  mpi_sgemm(m, n, k, leftMat, rightMat, resMat, rank, worldSize);
  
  gettimeofday(&stop, NULL);

  if (rank == 0) {
    cout << "mpi matmul: "
         << (stop.tv_sec - start.tv_sec) * 1000.0 +
                (stop.tv_usec - start.tv_usec) / 1000.0
         << " ms" << endl;

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < k; j++)
        if (int(resMat[i * k + j]) != n) {
          cout << resMat[i * k + j] << "error\n";
          exit(-1);
        }
    }
  }

  MPI_Finalize();
}
