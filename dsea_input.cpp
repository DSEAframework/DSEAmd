// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>

using namespace std;

int32_t DS::thread_input (int32_t n_super_cycle, int32_t myID, int32_t nProcs) {
	int32_t i_store=0;
	int64_t ta=MyGetTime();
	int64_t tb=-1;

	int64_t n_bytes_loaded=0;
	int64_t n_mol_loaded=0;
	int64_t current_nm_sum=-1;

	double perf_hist_10[10];
	for (int32_t i=0;i<10;i++) perf_hist_10[i]=-1;
	double perf_hist_100[100];
	for (int32_t i=0;i<100;i++) perf_hist_100[i]=-1;
	double perf_hist_1000[1000];
	for (int32_t i=0;i<1000;i++) perf_hist_1000[i]=-1;

	int32_t i_perf_hist_10=0;
	int32_t i_perf_hist_100=0;
	int32_t i_perf_hist_1000=0;

	for (int32_t i_supercycle=0;i_supercycle<n_super_cycle;i_supercycle++) {
		for (int32_t i_part=0;i_part<n_part;i_part++) {
			// get i_part to i_store
			while (stat_mem_in[i_store].state!=mem_state_free) {
				if (stat_mem_in[i_store].i_event!=-1) {
					// cout << "IN_event_" << stat_mem_in[i_store].i_event << endl;
					cudaError_t ces=cudaEventSynchronize (worker_event[stat_mem_in[i_store].i_event]);
					if (ces==cudaSuccess) {
						// cout << "IN:ready!" << endl;
						stat_mem_in[i_store].i_event=-1;
						stat_mem_in[i_store].i_part=-1;
						stat_mem_in[i_store].i_cycle=-1;
						stat_mem_in[i_store].n_use=0;
						stat_mem_in[i_store].state=mem_state_free;
					}
				}
			}	// wait till slot is free

			// cout << "IN_wait_complete" << endl;

			if ((i_supercycle==0)&&(myID==0)) {
				stat_mem_in[i_store].i_part=i_part;
				stat_mem_in[i_store].i_cycle=i_supercycle;

				int32_t my_tag=i_part*1000+968;

				MPI_Status status;
				int32_t ret=MPI_Recv((void *)d_in[i_store],store_size/sizeof(int32_t),MPI_INT,MPI_ANY_SOURCE,my_tag,MPI_COMM_WORLD,&status);

				if (ret==-1) {
					cout << "error_reading_file" << endl;
					// the program is in a broken state now... exit (todo)
				}
				else {
					// n_bytes_loaded+=dat_size[0];
					// dummy_sum is computed so that the compiler can not change the sequence
					// of writing data (1) and p_comm_in (2). the value written in (2) depends
					// the result of dummy_sum
					int32_t dummy_sum=ret;
					if (dummy_sum==2) {
						// stat_mem_in[i_store].i_part=i_part;	// (1)
						stat_mem_in[i_store].state=mem_state_ready;
					}
					else {
						// stat_mem_in[i_store].i_part=i_part;	//(2)
						stat_mem_in[i_store].state=mem_state_ready_b;
					}
//						pre_nm_sum+=BlockGetNM((char*)mycu.store_p_in[i_store]);
// 						current_nm_sum+=BlockGetNM((char*)mycu.store_p_in[i_store]);
// // 						cout << "L"<<i_part<< "_"<< i_store << "_" << dat_size[0] << endl;
				}
			}
			else {
					// receive part
				// cout << "MPI_Rec" << endl;
				int32_t source=myID-1;
				if (myID==0) source=nProcs-1;	// first rank receive from last rank

				int32_t tag=i_part;
				MPI_Status status;
				int32_t mpi_ret=MPI_Recv((void *)d_in[i_store],store_size/sizeof(int32_t),MPI_INT,source,tag,MPI_COMM_WORLD,&status);

				if (i_part==0) {
					tb=MyGetTime();
					if ((ta>0)&&(current_nm_sum>0)) {

						double seconds=(double)(tb-ta)/1000000.0;
						int32_t my_worker_sum=nProcs*n_worker;
						double rate=(double)current_nm_sum*(double)my_worker_sum/seconds;

						// record convergence data
						if (myID==0) {
							perf_hist_10[i_perf_hist_10]=rate;
							i_perf_hist_10++;
							i_perf_hist_10%=10;
							double perf_10=0;
							int n_perf_10=0;
							for(int32_t i=0;i<10;i++) {
								if (perf_hist_10[i]>0) {
									perf_10+=perf_hist_10[i];
									n_perf_10++;
								}
							}
							perf_10/=n_perf_10;

							perf_hist_100[i_perf_hist_100]=rate;
							i_perf_hist_100++;
							i_perf_hist_100%=100;
							double perf_100=0;
							int n_perf_100=0;
							for(int32_t i=0;i<100;i++) {
								if (perf_hist_100[i]>0) {
									perf_100+=perf_hist_100[i];
									n_perf_100++;
								}
							}
							perf_100/=n_perf_100;

							perf_hist_1000[i_perf_hist_1000]=rate;
							i_perf_hist_1000++;
							i_perf_hist_1000%=1000;
							double perf_1000=0;
							int n_perf_1000=0;
							for(int32_t i=0;i<1000;i++) {
								if (perf_hist_1000[i]>0) {
									perf_1000+=perf_hist_100[i];
									n_perf_1000++;
								}
							}
							perf_1000/=n_perf_1000;

							cout << rate << " " << perf_10 << " (10) " << perf_100 << " (100) " << perf_1000 << " (1000) " << " mol/s " << current_nm_sum << " nm " << seconds << " s " << endl;
							// cout << rate << " mol/s " << current_nm_sum << " nm " << seconds << " s " << (double)n_bytes_loaded/1.0e9/seconds << "GB/s #current_nm_sum" << endl;
						}
					}
					current_nm_sum=0;
					ta=tb;
				}

				if (mpi_ret==MPI_SUCCESS) {
					current_nm_sum+=get_block_nm_device((double*)d_in[i_store]);
					stat_mem_in[i_store].i_part=i_part;
					stat_mem_in[i_store].i_cycle=i_supercycle;
					stat_mem_in[i_store].state=mem_state_ready;
				}
			}

			i_store++;
			i_store%=n_store_in;
		}
		if ((i_supercycle==0)&&(myID==0)) {
			tb=MyGetTime();
			double seconds=(double)(tb-ta)/1000000.0;
			ta=tb;

			cout << "n_mol_loaded_" << n_mol_loaded << " bytes:" << n_bytes_loaded << " t: " << seconds << " " << (double)n_bytes_loaded/1.0e9/seconds << " GB/s" << endl;
		}
	}
	return 0;
}
