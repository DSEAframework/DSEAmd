// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>

using namespace std;

int32_t DS::thread_store_input (int32_t n_super_cycle, int32_t myID, int32_t nProcs) {
	int32_t i_store=0;

	for (int32_t i_supercycle=0;i_supercycle<n_super_cycle;i_supercycle++) {
		for (int32_t i_part=0;i_part<n_part;i_part++) {
			// get i_part to i_store
			while (stat_mem_store[i_store].state!=mem_state_free) {
				if (stat_mem_store[i_store].i_event!=-1) {
					// cout << "IN_event_" << stat_mem_store[i_store].i_event << endl;
					cudaError_t ces=cudaEventSynchronize (worker_event[stat_mem_store[i_store].i_event]);
					if (ces==cudaSuccess) {
						// cout << "IN:ready!" << endl;
						stat_mem_store[i_store].state=mem_state_free;
					}
				}
			}	// wait till slot is free

			// receive part
			// cout << "MPI_Rec" << endl;
			int32_t source=myID-1;
			if (myID==0) source=nProcs-1;	// first rank receive from last rank

			int32_t tag=i_part;
			MPI_Status status;
			int32_t mpi_ret=MPI_Recv((void *)h_store[i_store],store_size/sizeof(int32_t),MPI_INT,source,tag,MPI_COMM_WORLD,&status);

			if (mpi_ret==MPI_SUCCESS) {
				stat_mem_store[i_store].state=mem_state_ready;
			}

			i_store++;
			i_store%=n_store_host;
		}
	}
	return 0;
}
