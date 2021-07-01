
#include <mpi.h>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <iterator>
#include <limits>



inline int block_decompose(const int n, const int p, const int rank)
{
    return n / p + ((rank < n % p) ? 1 : 0);
}

inline int block_decompose(const int n, MPI_Comm comm)
{
    int rank, p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);
    return block_decompose(n, p, rank);
}

inline int block_decompose_by_dim(const int n, MPI_Comm comm, int dim)
{
    // get dimensions
    int dims[2];
    int periods[2];
    int coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);
    return block_decompose(n, dims[dim], coords[dim]);
}

template<typename T>
std::vector<T> read_binary_file(const char* filename)
{
    // get file size
    std::ifstream in(filename, std::ios::binary | std::ios::ate);
    if (!(in.good() && in.is_open()))
        throw std::runtime_error(std::string("Couldn't open file ") + filename);
    std::size_t nbytes = in.tellg();
    in.close();
    // open again, this time at the beginning
    std::ifstream infile(filename, std::ios::binary);

    // create vector of the correct size
    std::size_t n = nbytes / sizeof(T);
    std::vector<T> result(n);
    // read file into vector:
    infile.read((char*) &result[0], nbytes);
    return result;
}


template<typename T>
void write_binary_file(const char* filename, const std::vector<T>& data)
{
    // open file in binary mode
    std::ofstream out(filename, std::ios::binary);
    if (!(out.good() && out.is_open()))
        throw std::runtime_error(std::string("Couldn't open file ") + filename);
    // write data to file
    out.write(reinterpret_cast<const char*>(&data[0]), sizeof(T)*data.size());
    out.close();
}


// add code begin
/*
inline void GaussianElimination(){
    
}
*/
void jacobi_Seq(std::vector<double> &A,std::vector<double> &b ,std::vector<double> &x,int n){
    const double  epsilon = 1e-4;
    //Ax = b   (n*n)*(n*1) == n*1
    for(int i=0;i<n;i++){
        x[i] = b[i]/A[i*n+i];//initial estimation
    }
    double  diff = 0;
    std::vector<double> newx(b);
    do{
        diff = 0;
        for(int i=0;i<n;i++){
            
            newx[i] = x[i];
            for(int j=0;j<n;j++){
                if(i!=j)
                    newx[i] = newx[i]-A[i*n+j]*x[i];
            }
            newx[i] = newx[i]/A[i*n+n];            
        }
        for(int i=0;i<n;i++){
            diff = std::max(diff,abs(x[i]-newx[i]));
            x[i] = newx[i];
        }
    }while(diff > epsilon);

    printf("diff:%lf\n",diff);
}
// add code end



int main(int argc, char *argv[])
{
   // set up MPI
   MPI_Init(&argc, &argv);

   // get communicator size
   MPI_Comm comm = MPI_COMM_WORLD;
   int p;
   MPI_Comm_size(comm, &p);

   // Ax = b
   std::vector<double> A;
   std::vector<double> b;
   std::vector<double> x;



   std::string outfile_name;

   std::string fileA;
   std::string fileB;

    // get output filename
  outfile_name = std::string(argv[3]);

  int n;

   // start timer
   //   we omit the file loading and argument parsing from the runtime
   //   timings, we measure the time needed by the processor with rank 0
   struct timespec t_start, t_end;

  // get the dimensions
  int q = (int)sqrt(p);
  if (p != q*q)
  {
     throw std::runtime_error("The number of MPI processes must be a perfect square");
  }
  // create 2D cartesian grid for the processors (enable reordering)
  MPI_Comm grid_comm;
  int dims[2] = {q, q};
  int periods[2] = {0, 0};
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);
  // get the rank of process with coordinates (0,0)
  int rank00;
 
  int myrank;
  int coords[2] = {0, 0};
  MPI_Cart_rank(grid_comm, coords, &rank00);
  MPI_Comm_rank(grid_comm, &myrank);

  // process (0,0) loads input
  if (myrank == rank00)
  {

    // read input from file
    A = read_binary_file<double>(argv[1]);
    b = read_binary_file<double>(argv[2]);
    n = b.size();
    if (A.size() != n*(size_t)n)
    {
       throw std::runtime_error("The input dimensions are not matching");
    }
        
  }
  MPI_Bcast(&n, 1, MPI_INT, rank00, grid_comm);

  // start timer
  clock_gettime(CLOCK_MONOTONIC,  &t_start);


  if (myrank == rank00)
     x = std::vector<double>(n);


// add code begin
    //jacobi_Seq(A,b,x,n);
    //deal with this problem with p processors
    const double eps = 1e-4;
    const int MAX_ITERATIONS = 1e4;
    
    if(myrank == rank00){
        int sended_rows=block_decompose(n,p,myrank);
        for(int i=1;i<p;i++){
            int k = block_decompose(n,p,i);
            MPI_Send(&A[sended_rows*n],k*n,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
            MPI_Send(&b[sended_rows],n,MPI_DOUBLE,i,1,MPI_COMM_WORLD);
        }
    }else{
            int k = block_decompose(n,p,myrank);
            A.resize(k*n);
            b.resize(n);
            x.resize(n);
            MPI_Recv(&A[0],k*n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(&b[0],n,MPI_DOUBLE,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    //确保所有进程数据分配完毕
    MPI_Barrier(MPI_COMM_WORLD);
    //processors from 0 to p-1 
    int iter = 0;
    int k = block_decompose(n,p,myrank);
    double diff = 0;
    std::vector<double> newx(b);
    
    int *recvcounts  = new int[p];
    int *displs = new int[p];
    recvcounts[0] = block_decompose(n,p,0);
    displs[0] = 0;
    for(int i=1;i<p;i++){
        recvcounts[i] = block_decompose(n,p,i);
        displs[i] = recvcounts[i-1]+displs[i-1];
    }
    do{
        //x.assign(newx.begin(),newx.end());
        for(int i=0;i<k;i++){
            double sigma = 0;
            for(int j=0;j<n;j++){
                if(i!=j)
                    sigma += A[i*n+j]*x[j];
            }
            newx[i] = (b[i]- sigma)/A[i*n+i];
        }
        //将newx的数据组合同步到所有进程
        // MPI_Allgatherv(const void * sendbuf，int sendcount，MPI_Datatype sendtype，
        //                 void * recvbuf，const int * recvcounts，const int * displs，
        //                 MPI_Datatype recvtype，MPI_Comm comm)
        //搞定这一步就万事大吉了哈哈哈哈        
        MPI_Allgatherv(&newx[displs[myrank]],recvcounts[myrank],MPI_DOUBLE,&x[0],recvcounts,displs,MPI_DOUBLE,MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);
        diff = 0;
        for(int i=0;i<n;i++){
            diff = std::max(diff,abs(x[i]-newx[i]));
        }
        iter++;
        
    }while(iter < MAX_ITERATIONS && diff > eps);
// add code end

  if (myrank == rank00)
  {
     // get time
     clock_gettime(CLOCK_MONOTONIC,  &t_end);
     // time in seconds
     double time_secs = (t_end.tv_sec - t_start.tv_sec)
        + (double) (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
     // output time
     std::cerr << time_secs << std::endl;
     // write output

        write_binary_file(outfile_name.c_str(), x);
     
  }
   
   // finalize MPI
   MPI_Finalize();
   return 0;
}
