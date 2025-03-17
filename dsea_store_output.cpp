// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>

using namespace std;

int32_t DS::thread_store_output (int32_t n_super_cycle, int32_t myID, int32_t nProcs) {
	int32_t i_store=0;

	for (int32_t i_supercycle=0;i_supercycle<n_super_cycle;i_supercycle++) {

		// regular cycle
		for (int32_t i_part=0;i_part<n_part;i_part++) {
			// cout << "OUT:start!_" << i_store << "_" << stat_mem_out[i_store].i_event << endl;
			while (stat_mem_store[i_store].i_event==mem_state_ready) {}

			int32_t dest=myID+1;
			if (dest==nProcs) dest=0;
			int32_t tag=i_part;
			// cout << "MPI_Send_" << tag << "_" << dest << endl;

			int32_t send_size=get_block_size_device((double*)h_store[i_store]);
			int32_t res=MPI_Send((const void*)h_store[i_store],send_size/sizeof(int32_t),MPI_INT,dest,tag,MPI_COMM_WORLD);

			if (res==MPI_SUCCESS) {
				stat_mem_store[i_store].state=mem_state_free;
			}

			i_store++;
			i_store%=n_store_host;
		}
	}

	return 0;
}
