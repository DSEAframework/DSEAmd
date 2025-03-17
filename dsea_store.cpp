// Data Streaming for Explicit Algorithms - DSEA

#include <iostream>
#include <mpi.h>
#include <dsea.h>
using namespace std;

int32_t main(int argc, char ** argv) {

	if (argc!=2) {
		cout << "usage: dsea-sore n_cycles" << endl;
		return -1;
	}
	int32_t n_super_cycle=atoi(argv[1]);
	// int32_t n_host_store=atoi(argv[2]);

	int32_t tmp_rank=0;		// in case MPI is not used
	int32_t tmp_nProcs=1;	// in case MPI is not used
	// char ProcessorName [1000];

	int32_t provided=-1;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
	MPI_Comm_rank(MPI_COMM_WORLD,&tmp_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&tmp_nProcs);

	int32_t myID=tmp_rank;
	int32_t nProcs=tmp_nProcs;

	// for (int32_t i=0;i<nProcs;i++) {
 		// if (myID==i) cout << "INFO: rank " << myID << " running on: " << endl;// ProcessorName << /*" " << myNUMAnode <<      " " << myIDhost << " " << myIdNUMA << " " << HostMaster << " " << NUMANodeMaster << " " << myJob <<*/ endl;
 			// MPI_Barrier(MPI_COMM_WORLD);
 	// }

	// for (int32_t i_super_cycle=0;i_super_cycle<n_super_cycle;i_super_cycle++) {
	// 	for (int32_t i_part=0;i_part<my_n_part;i_part++) {
			
	// 	}
	// }

	DS ds(0,0,my_n_part,0,0,myID,nProcs);

#pragma omp parallel default (none) num_threads(2) shared (cout) \
shared (ds,n_super_cycle) \
shared (myID,nProcs)
{
	#pragma omp single nowait	/* input thread */
	{
		ds.thread_store_input(n_super_cycle,myID,nProcs);
	}

	#pragma omp single nowait	/* output thread */
	{
		ds.thread_store_output(n_super_cycle,myID,nProcs);
	}
}
	MPI_Finalize();
}
