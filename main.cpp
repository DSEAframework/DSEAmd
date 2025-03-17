// Data Streaming for Explicit Algorithms - DSEA

#include <iostream>
#include <mpi.h>
#include <dsea.h>
using namespace std;

// #define my_n_worker			1

int32_t main(int argc, char ** argv) {

	if (argc!=3) {
		cout << "usage: dsea n_worker n_cycles n_rails" << endl;
//		return -1;
	}
	int32_t my_n_worker=atoi(argv[1]);
	int32_t n_super_cycle=atoi(argv[2]);
	int32_t my_n_rails=atoi(argv[3]);

	// cout << "order: " << order_in << " " << order_out << endl;
	// cout << "n_worker: " << my_n_worker << endl;
	// cout << "n_part: " << my_n_part << endl;
	// cout << endl;

	int32_t tmp_rank=0;		// in case MPI is not used
	int32_t tmp_nProcs=1;	// in case MPI is not used
	char ProcessorName [1000];

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

	int32_t igpu=0;//myID % 8;
	int32_t my_order_in=1;
	int32_t my_order_out=1;

	DS ds(igpu,my_n_worker,my_n_part,my_order_in,my_order_out,myID,nProcs,my_n_rails);

#pragma omp parallel default (none) num_threads(4) shared (cout) \
shared (ds,n_super_cycle) \
shared (myID,nProcs,argc,argv)
{
	ds.CudaDummy();

	#pragma omp master			/* worker thread */
	{
		// cout << "thread_master_start_" << myID << endl;
		ds.thread_main(n_super_cycle,1,1,myID);
		// cout << "thread_master_done_" << myID << endl;
	}

	#pragma omp single nowait	/* input thread */
	{
		// cout << "comm_thread_IN_start" << myID << endl;
#ifndef MRUCX_REC
		// mpi only version
		// cout << "MRUCX_OFF" << endl;
		ds.thread_input(n_super_cycle,myID,nProcs);
#else 
		// cout << "MRUCX_REC" << endl;
		ds.thread_input_ucx(argc,argv,n_super_cycle,myID,nProcs);
#endif
		// cout << "comm_thread_IN_done" << myID << endl;
	}

	#pragma omp single nowait	/* output thread */
	{
		// cout << "comm_thread_OUT_start" << myID << endl;
#ifndef MRUCX_SEND
		// cout << "MRUCX_OFF" << endl;
		ds.thread_output(n_super_cycle,myID,nProcs);
#else
		// cout << "MRUCX_SEND" << endl;
		ds.thread_output_ucx(argc,argv,n_super_cycle,myID,nProcs);
#endif
		// cout << "comm_thread_OUT_done" << myID << endl;
	}

	#pragma omp single nowait	/* storage thread */
	{
		// cout << "comm_thread_OUT_start" << myID << endl;
		ds.thread_storage(n_super_cycle,myID,nProcs);
		// cout << "comm_thread_OUT_done" << myID << endl;
	}

}

	MPI_Finalize();
}
