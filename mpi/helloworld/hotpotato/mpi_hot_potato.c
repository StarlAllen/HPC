
#include <stdio.h>
#include <mpi.h>

#define first 0

int main(int argc,char ** argv){
	int numprocs,world_rank,world_size,namelen;
	char potato[] = "Greetings";
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;

	MPI_Init(&argc,&argv);
	//MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	MPI_Get_processor_name(processor_name,&namelen);
        if(world_rank==first){
	//               buff   item item_size  src/dest  tag world
		MPI_Send(potato,1,  MPI_INT,    1,        0, MPI_COMM_WORLD);
              	printf("Greetings!\n");
		//MPI_Recv(&potato,1,MPI_INT,world_size-1,0,MPI_COMM_WORLD,&status);
	
	}else if(world_rank < world_size-1){
               MPI_Recv(potato,1,MPI_INT,world_rank-1,0,MPI_COMM_WORLD,&status);
	       MPI_Send(potato,1,MPI_INT,(world_rank+1)%world_size,0,MPI_COMM_WORLD);
       	       printf("Message send over from processor %s,rank %d out of %d processors\n",
			       processor_name,world_rank,world_size);
	}else{
		MPI_Recv(&potato,1,MPI_INT,world_rank-1,0,MPI_COMM_WORLD,&status);
		printf("Message send over from processor %s,rank %d out of %d processors\n",
				processor_name,world_rank,world_size);
	}

	MPI_Finalize();
	return 0;
}
