// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include "src/ucx_multirail.h"

#include <chrono>	// sleep
#include <thread>	// sleep

using namespace std;

int32_t DS::thread_input_ucx (int argc, char ** argv,int32_t n_super_cycle, int32_t myID, int32_t nProcs) {
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

	// init ucx dual rail
    cudaSetDevice(0);

    ucs_status_t status;
    ucx_mr_context_t mr_ctx;
    mr_ctx.server_addr = NULL;

    parse_opts(&mr_ctx, argc, argv);

    status = ucx_mr_setup(&mr_ctx);
    if (status != UCS_OK) {
        printf("There was a problem!\n");
    }
    ucx_mr_test_connection(&mr_ctx);

	std::this_thread::sleep_for(std::chrono::milliseconds(2000));

	cudaSetDevice(1);
	CUdeviceptr tmp_store1;
	cudaMalloc((void**)&tmp_store1,size_device);
	cudaSetDevice(2);
	CUdeviceptr tmp_store2;
	cudaMalloc((void**)&tmp_store2,size_device);
	cudaSetDevice(3);
	CUdeviceptr tmp_store3;
	cudaMalloc((void**)&tmp_store3,size_device);
	cudaSetDevice(0);
	CUdeviceptr tmp_store0;
	cudaMalloc((void**)&tmp_store0,size_device);

	// {
	// // cout << "thread_input_ucx:init_done" << endl;

	// int64_t t_single_a=MyGetTime();
	int n_loop=20;
	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x133;
	// 	ucx_mr_single_recv(&mr_ctx, 0, tag, (void*)tmp_store0, size_device, UCS_MEMORY_TYPE_CUDA, 0);
	// }
	// int64_t t_single_b=MyGetTime();
	// double t_single=(t_single_b-t_single_a)/(double)n_loop/1e6;
	// double rate_single=(double)size_device/t_single;
	// rate_single/=1e9;
	// cout << "rate_single_rail_0:_" << rate_single << endl;
	// // cudaSetDevice(0);



	// t_single_a=MyGetTime();
	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x233;
	// 	ucx_mr_single_recv(&mr_ctx, 1, tag, (void*)tmp_store1, size_device, UCS_MEMORY_TYPE_CUDA, 1);
	// }
	// t_single_b=MyGetTime();
	// t_single=(t_single_b-t_single_a)/(double)n_loop/1e6;
	// rate_single=(double)size_device/t_single;
	// rate_single/=1e9;
	// cout << "rate_single_rail_1:_" << rate_single << endl;
	// // cudaSetDevice(0);


	// t_single_a=MyGetTime();
	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x333;
	// 	ucx_mr_single_recv(&mr_ctx, 2, tag, (void*)tmp_store2, size_device, UCS_MEMORY_TYPE_CUDA, 2);
	// }
	// t_single_b=MyGetTime();
	// t_single=(t_single_b-t_single_a)/(double)n_loop/1e6;
	// rate_single=(double)size_device/t_single;
	// rate_single/=1e9;
	// cout << "rate_single_rail_2:_" << rate_single << endl;
	// // cudaSetDevice(0);


	// t_single_a=MyGetTime();
	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x433;
	// 	ucx_mr_single_recv(&mr_ctx, 3, tag, (void*)tmp_store3, size_device, UCS_MEMORY_TYPE_CUDA, 3);
	// }
	// t_single_b=MyGetTime();
	// t_single=(t_single_b-t_single_a)/(double)n_loop/1e6;
	// rate_single=(double)size_device/t_single;
	// rate_single/=1e9;
	// cout << "rate_single_rail_3:_" << rate_single << endl;
	// // cudaSetDevice(0);




	// // // cout << "thread_input_ucx:single rec done" << endl;
	// cout << "starting dual rail benchmark..." << endl;

	// int64_t t_dual_a=MyGetTime();
	// n_loop=1000;
	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x533;
	// 	int element_size = sizeof(int32_t);
	// 	float split_ratio=0.5;
	// 	int32_t ucx_ret=ucx_mr_split_recv(&mr_ctx, tag, split_ratio, element_size,(void*)d_in[0], size_device, UCS_MEMORY_TYPE_CUDA, 0,(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1);
	// 	// int32_t ucx_ret=ucx_mr_dual_split_recv(&mr_ctx, tag, split_ratio, element_size,(void*)d_in[i_store], size_device, UCS_MEMORY_TYPE_CUDA, 0,(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,1);
	// }
	// int64_t t_dual_b=MyGetTime();
	// double t_dual=(t_dual_b-t_dual_a)/(double)n_loop/1e6;
	// double rate_dual=(double)size_device/t_dual;
	// rate_dual/=1e9;
	// cout << "rate_dual:_" << rate_dual << endl;


	int64_t t_triple_a=MyGetTime();
	n_loop=200;
	for (int i=0;i<n_loop;i++) {
		// cout << i << endl;
		ucp_tag_t tag = 0x733;
		int element_size = sizeof(int32_t);
		float split_ratio=0.66;
		
		int32_t ucx_ret=ucx_mr_tripple_split_recv_simple(&mr_ctx, tag, split_ratio, element_size,(void*)d_in[0], size_device, UCS_MEMORY_TYPE_CUDA, 0,
		(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,
		(void*)tmp_store2, UCS_MEMORY_TYPE_CUDA, 2);
		if (ucx_ret!=UCS_OK) {
			cout << "prob_rec_a3" << endl;
		}
	}
	int64_t t_tripple_b=MyGetTime();
	double t_tripple=(t_tripple_b-t_triple_a)/(double)n_loop/1e6;
	double rate_tripple=(double)size_device/t_tripple;
	rate_tripple/=1e9;
	cout << "rate_tripple:_" << rate_tripple << endl;



	// int64_t t_quad_a=MyGetTime();
	// n_loop=2000;
	// for (int i=0;i<n_loop;i++) {

	// 	if (i%100==0) cout << i << endl;
	// 	ucp_tag_t tag = 0x633;
	// 	int element_size = sizeof(int32_t);
	// 	float split_ratio=0.75;
		
	// 	int32_t ucx_ret=ucx_mr_quad_split_recv_simple(&mr_ctx, tag, split_ratio, element_size,(void*)d_in[0], size_device, UCS_MEMORY_TYPE_CUDA, 0,
	// 	(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,
	// 	(void*)tmp_store2, UCS_MEMORY_TYPE_CUDA, 2,
	// 	(void*)tmp_store3, UCS_MEMORY_TYPE_CUDA, 3);
	// 	if (ucx_ret!=UCS_OK) {
	// 		cout << "prob_rec_a4" << endl;
	// 	}
	// }
	// int64_t t_quad_b=MyGetTime();
	// double t_quad=(t_quad_b-t_quad_a)/(double)n_loop/1e6;
	// double rate_quad=(double)size_device/t_quad;
	// rate_quad/=1e9;
	// cout << "rate_quad:_" << rate_quad << endl;



    // cout << "thread_input_ucx:ready_to_go" << endl;
	// }

	// return 0;

	cudaSetDevice(0);


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
				// cout << "wait_" << i_part << endl;
				int32_t ret=MPI_Recv((void *)d_in[i_store],store_size/sizeof(int32_t),MPI_INT,MPI_ANY_SOURCE,my_tag,MPI_COMM_WORLD,&status);
				// cout << "got_" << i_part << endl;

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

				// int32_t tag=i_part;
				// MPI_Status status;
				// cout << "debug_MPI_Rec_pre_" << tag << endl;

				// int32_t mpi_ret=MPI_Recv((void *)d_in[i_store],store_size/sizeof(int32_t),MPI_INT,source,tag,MPI_COMM_WORLD,&status);
				// cout << "debug_MPI_Rec_post_" << tag << endl;

				int element_size = sizeof(int32_t);
				// cout << "pre_ucx_mr_split_recv" << endl;

				int myrails=n_rails;
				int32_t ucx_ret=-1;
				if ((i_supercycle==1)&&(myID==0)&&(i_part==0)) cout << "rails_used:" << myrails << endl;

				if (myrails==1) {
					// 1 rail
					ucp_tag_t tag = 0x51;
					ucx_ret=ucx_mr_single_recv(&mr_ctx, 0, tag, (void*)d_in[i_store], store_size, UCS_MEMORY_TYPE_CUDA, 0);
				}
				else if (myrails==2) {
					// 2 rails
					ucp_tag_t tag = 0x52;
					float split_ratio=0.5;
					ucx_ret=ucx_mr_split_recv(&mr_ctx, tag, split_ratio, element_size,(void*)d_in[i_store], store_size, UCS_MEMORY_TYPE_CUDA, 0,
					(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1);
					// ucx_ret=ucx_mr_dual_split_recv(&mr_ctx, tag, split_ratio, element_size,(void*)d_in[i_store], store_size, UCS_MEMORY_TYPE_CUDA, 0,(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,4);
				}
				else if (myrails==3) {
					ucp_tag_t tag = 0x53;
					float split_ratio=0.66;
					ucx_ret=ucx_mr_tripple_split_recv_simple(&mr_ctx, tag, split_ratio, element_size,(void*)d_in[i_store], store_size, UCS_MEMORY_TYPE_CUDA, 0,
					(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,
					(void*)tmp_store2, UCS_MEMORY_TYPE_CUDA, 2);
					// cout << "rec_check" << endl;
					// block_check_device((double*)d_in[i_store],123);
					// cout << "rec_check_done" << endl;

				}
				else if (myrails==4) {
					ucp_tag_t tag = 0x54;
					float split_ratio=0.75;
					ucx_ret=ucx_mr_quad_split_recv_simple(&mr_ctx, tag, split_ratio, element_size,(void*)d_in[i_store], store_size, UCS_MEMORY_TYPE_CUDA, 0,
					(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,
					(void*)tmp_store2, UCS_MEMORY_TYPE_CUDA, 2,
					(void*)tmp_store3, UCS_MEMORY_TYPE_CUDA, 3);
					// cout << "rec_check" << endl;
					// block_check_device((double*)d_in[i_store],123);
					// cout << "rec_check_done" << endl;

				}



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
						}
					}
					current_nm_sum=0;
					ta=tb;
				}

				// if (mpi_ret==MPI_SUCCESS) {
				if (ucx_ret==0) {
					current_nm_sum+=get_block_nm_device((double*)d_in[i_store]);
					stat_mem_in[i_store].i_part=i_part;
					stat_mem_in[i_store].i_cycle=i_supercycle;
					stat_mem_in[i_store].state=mem_state_ready;
				}
				else {
					cout << "fail:_ucx_mr_split_recv_" << ucx_ret << endl;
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

    ucx_mr_cleanup(&mr_ctx, FULL);

	return 0;
}
