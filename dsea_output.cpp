// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include <fstream>


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


int32_t DS::thread_output (int32_t n_super_cycle, int32_t myID, int32_t nProcs) {
	int32_t i_store=0;

	int64_t n_mol_stored=0;



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
					int32_t tag=i_part;
					// cout << "MPI_Send_" << tag << "_" << dest << endl;

					int32_t send_size=get_block_size_device((double*)d_out[i_store]);

					int32_t res=MPI_Send((const void*)d_out[i_store],send_size/sizeof(int32_t),MPI_INT,dest,tag,MPI_COMM_WORLD);

					if (res==MPI_SUCCESS) {
						stat_mem_out[i_store].i_event=-1;
						stat_mem_out[i_store].i_part=-1;
						stat_mem_out[i_store].n_use=0;
						stat_mem_out[i_store].state=mem_state_free;
					}

					i_store++;
					if (i_store==n_store_out) i_store=0;
				}
			}
		}
	}

	return 0;
}
