// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include <fstream>
#include "src/ucx_multirail.h"

#include <chrono>	// sleep
#include <thread>	// sleep

using namespace std;

int64_t DS::MemToFile(int64_t * dat, int64_t size, char * FileName, int32_t newfile) {
        // debug
        std::string txt; 
        std::string debug_filename;

//      cout << dat << " " << n << " " << newfile << endl;
     
        ofstream ofs; 
        if (newfile==1) {    
                ofs.open(FileName, ios::out | ios::binary);
//              cout << FileName << " new " << n << endl;
        }    
        else if (newfile==0) {    
                ofs.open(FileName, ios::out | ios::binary | ios::app);
//              cout << FileName << " app " << n << endl;
        }    
        else {    
                cout << "invalid MemToFile" << endl;
                return -1;
        }    
     
        if (ofs) {    
                ofs.write((char*)dat,size);
                if (ofs) {    
//                      cout << n << " bytes written" << endl;
                        ofs.close();
                        return size;
                }    
                else {    
                        cout << "problem writing file: " << FileName << " " << size << endl;
                        cout << "aborting..." << endl;
//                      exit(-1);
                }    
                ofs.close();
        }    
        else {    
                cout << "problem opening file: " << FileName << endl;
        }    
        return -1;
}


int32_t DS::thread_output_ucx (int argc, char ** argv, int32_t n_super_cycle, int32_t myID, int32_t nProcs) {
	int32_t i_store=0;
	int64_t n_mol_stored=0;

	std::this_thread::sleep_for(std::chrono::milliseconds(4000));

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
	// cout << "thread_output_ucx:init_done" << endl;

	int n_loop=20;
	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x133;
	// 	ucx_mr_single_send(&mr_ctx, 0, tag,(void*)tmp_store0, size_device, UCS_MEMORY_TYPE_CUDA, 0);
	// }
	// // cudaSetDevice(0);
	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x233;
	// 	ucx_mr_single_send(&mr_ctx, 1, tag,(void*)tmp_store1, size_device, UCS_MEMORY_TYPE_CUDA, 1);
	// }
	// // cudaSetDevice(0);

	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x333;
	// 	ucx_mr_single_send(&mr_ctx, 2, tag,(void*)tmp_store2, size_device, UCS_MEMORY_TYPE_CUDA, 2);
	// }
	// // cudaSetDevice(0);

	// for (int i=0;i<n_loop;i++) {
	// 	ucp_tag_t tag = 0x433;
	// 	ucx_mr_single_send(&mr_ctx, 3, tag,(void*)tmp_store3, size_device, UCS_MEMORY_TYPE_CUDA, 3);
	// }
	// // cudaSetDevice(0);

	// cout << "thread_output_ucx:single send done" << endl;



	// n_loop=1000;
	// for (int i=0;i<n_loop;i++) {
	// 	int element_size = sizeof(int32_t);
	// 	float split_ratio=0.5;
	// 	ucp_tag_t tag = 0x533;
	// 	int32_t ucx_ret=ucx_mr_split_send(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], size_device, UCS_MEMORY_TYPE_CUDA, 0,(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1);
	// 	// int32_t ucx_ret=ucx_mr_dual_split_send(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], size_device, UCS_MEMORY_TYPE_CUDA, 0,(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,1);
	// }

	n_loop=200;


	for (int i=0;i<n_loop;i++) {
		int element_size = sizeof(int32_t);
		float split_ratio=0.66;
		ucp_tag_t tag = 0x733;
		int32_t ucx_ret=ucx_mr_tripple_split_send_simple(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], size_device, UCS_MEMORY_TYPE_CUDA, 0,
						(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,
						(void*)tmp_store2, UCS_MEMORY_TYPE_CUDA, 2);
		if (ucx_ret!=UCS_OK) {
			cout << "prob_send_a3_" << ucx_ret << endl;
		}
	}


	// for (int i=0;i<n_loop;i++) {
	// 	int element_size = sizeof(int32_t);
	// 	float split_ratio=0.75;
	// 	ucp_tag_t tag = 0x633;
	// 	int32_t ucx_ret=ucx_mr_quad_split_send_simple(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], size_device, UCS_MEMORY_TYPE_CUDA, 0,
	// 					(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,
	// 					(void*)tmp_store2, UCS_MEMORY_TYPE_CUDA, 2,
	// 					(void*)tmp_store3, UCS_MEMORY_TYPE_CUDA, 3);
	// 	if (ucx_ret!=UCS_OK) {
	// 		cout << "prob_send_a4" << endl;
	// 	}
	// 	// int32_t ucx_ret=ucx_mr_dual_split_send(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], size_device, UCS_MEMORY_TYPE_CUDA, 0,(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,1);
	// }

	cout << "thread_output_ucx:ready_to_go" << endl;
	// }

// return 0;


	cudaSetDevice(0);


	for (int32_t i_supercycle=0;i_supercycle<n_super_cycle;i_supercycle++) {

		if ((i_supercycle==n_super_cycle-1)&&(myID==nProcs-1)) {
			// last cycle in last MPI rank stores output
			char * my_block = MyNewCatchChar(__FILE__,__LINE__,store_size);

			for (int32_t i_part=0;i_part<n_part;i_part++) {
				// cout << "OUT:start!_" << i_store << "_" << stat_mem_out[i_store].i_event << endl;
				while (stat_mem_out[i_store].i_event==-1) {}	// wait for event

				// cout << "OUT:wait!_" << i_store << "_" << stat_mem_out[i_store].i_event << endl;
				cudaError_t ces=cudaEventSynchronize (worker_event[stat_mem_out[i_store].i_event]);

				if (ces==cudaSuccess) {
					// download block from GPU
					cudaMemcpy((void*)my_block,(const void*)d_out[i_store],store_size,cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);
					// block_check((double*)my_block,2);
					// block_markcheck((double*)my_block);
					n_mol_stored+=block_get_nm((double*)my_block);

					if (false) {
						string FileName;
						// FileName="rho_0.75_T_0.8/";
						FileName="rho_0.5_T_1.5/";
#if defined case_30
						FileName+="case_30/stream_out_";
#elif defined case_63
						FileName+="case_63/stream_out_";
#elif defined case_136
						FileName+="case_136/stream_out_";
#elif defined case_293
						FileName+="case_293/stream_out_";
#elif defined case_630
						FileName+="case_630/stream_out_";
#elif defined case_900
						FileName+="ws_listcase_900/stream_out_";
#else
						FileName+="?/stream_out_";
#endif

						FileName+=to_string(i_part);

						int64_t block_size=get_block_size_host((double*)my_block);
						MemToFile(&block_size,sizeof(int64_t),(char*)FileName.c_str(),1);

						MemToFile((int64_t*)my_block,block_size,(char*)FileName.c_str(),0);

						int64_t check_int=CHECK_INT;
						MemToFile(&check_int,sizeof(int64_t),(char*)FileName.c_str(),0);
					}

					// cout << "OUT:ready!_" << i_store << endl;
					stat_mem_out[i_store].i_event=-1;
					stat_mem_out[i_store].i_part=-1;
					stat_mem_out[i_store].n_use=0;
					stat_mem_out[i_store].state=mem_state_free;

					i_store++;
					if (i_store==n_store_out) i_store=0;
				}
			}

			cout << "n_mol_stored_" << n_mol_stored << endl;

			
			cudaMemcpy((void*)my_block,(const void*)d_debug,size_debug,cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);
			for (int i=0;i<12;i++) cout << i << "_" << ((int32_t *)my_block)[i] << endl;

			delete [] my_block;

		}
		else {
			// regular cycle
			for (int32_t i_part=0;i_part<n_part;i_part++) {
				// cout << "OUT:start!_" << i_store << "_" << stat_mem_out[i_store].i_event << endl;
				while (stat_mem_out[i_store].i_event==-1) {}	// wait for event

				// cout << "OUT:wait!_" << i_store << "_" << stat_mem_out[i_store].i_event << endl;
				cudaError_t ces=cudaEventSynchronize (worker_event[stat_mem_out[i_store].i_event]);
				if (ces==cudaSuccess) {
					int32_t dest=myID+1;
					if (dest==nProcs) dest=0;
					// int32_t tag=i_part;
					// cout << "MPI_Send_" << tag << "_" << dest << endl;

					int32_t send_size=get_block_size_device((double*)d_out[i_store]);

// cout << "debug_MPI_Send_pre_" << tag << endl;
					// int32_t res=MPI_Send((const void*)d_out[i_store],send_size/sizeof(int32_t),MPI_INT,dest,tag,MPI_COMM_WORLD);
// cout << "debug_MPI_Send_post_" << tag << endl;


					// //   void **buffer = mr_bench_ctx->send_buffer;
					// // ucs_memory_type_t mem_type = mr_ctx.mem_type;
					// //   size_t length = mr_bench_ctx->msg_size;

					int element_size = sizeof(int32_t);

					int myrails=n_rails;
					int32_t ucx_ret=-1;

					// cout << "pre_ucx_mr_split_send" << endl;
					if (myrails==1) {
						// 1 rail
						ucp_tag_t tag = 0x51;
						ucx_ret=ucx_mr_single_send(&mr_ctx, 0, tag,(void*)d_out[i_store], send_size, UCS_MEMORY_TYPE_CUDA, 0);
					}
					else if (myrails==2) {
						// 2 rails
						ucp_tag_t tag = 0x52;
						float split_ratio=0.5;
						ucx_ret=ucx_mr_split_send(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], send_size, UCS_MEMORY_TYPE_CUDA, 0,
						(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1);
						// ucx_ret=ucx_mr_dual_split_send(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], send_size, UCS_MEMORY_TYPE_CUDA, 0,(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,4);
					}
					else if (myrails==3) {
						ucp_tag_t tag = 0x53;
						float split_ratio=0.66;
						// cout << "send_check" << endl;
						// block_check_device((double*)d_out[i_store],123);
						// cout << "send_check_done" << endl;

						ucx_ret=ucx_mr_tripple_split_send_simple(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], send_size, UCS_MEMORY_TYPE_CUDA, 0,
						(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,
						(void*)tmp_store2, UCS_MEMORY_TYPE_CUDA, 2);

					}
					else if (myrails==4) {
						ucp_tag_t tag = 0x54;
						float split_ratio=0.75;
						// cout << "send_check" << endl;
						// block_check_device((double*)d_out[i_store],123);
						// cout << "send_check_done" << endl;

						ucx_ret=ucx_mr_quad_split_send_simple(&mr_ctx, tag, split_ratio, element_size,(void*)d_out[i_store], send_size, UCS_MEMORY_TYPE_CUDA, 0,
						(void*)tmp_store1, UCS_MEMORY_TYPE_CUDA, 1,
						(void*)tmp_store2, UCS_MEMORY_TYPE_CUDA, 2,
						(void*)tmp_store3, UCS_MEMORY_TYPE_CUDA, 3);

					}


					// if (res==MPI_SUCCESS) {
					if (ucx_ret==0) {
						stat_mem_out[i_store].i_event=-1;
						stat_mem_out[i_store].i_part=-1;
						stat_mem_out[i_store].n_use=0;
						stat_mem_out[i_store].state=mem_state_free;
					}
					else {
						cout << "fail:_ucx_mr_split_send_" << ucx_ret << endl;
					}

					i_store++;
					if (i_store==n_store_out) i_store=0;
				}
			}
		}
	}

    ucx_mr_cleanup(&mr_ctx, FULL);

	return 0;
}
