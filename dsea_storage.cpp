// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include <string>
#include <fstream>

using namespace std;

int32_t DS::FileFieldsToMem(int64_t * dat, char * FileName, int32_t iMax, int64_t * size) {
	int64_t int8_JCHECK=CHECK_INT;
	
	ifstream ifs;
	bool check_int_ok=false;
	
	ifs.open(FileName, ios::in | ios::binary);
	if (!ifs) {
		cout << "problem reading from file" << endl;
		cout << FileName << endl;
		return (-1);
	}
	else {
// 		if (myID==0) cout << FileName << " " << sizeof(streamsize) << endl;
	}

	int32_t iPointer=0;
	int64_t nbytes=0;
	int64_t bytes_read=0;	// number of bytes read

	for (int32_t i=0;i<16;i++) {
		*(dat+i) = (int64_t) 0;
		*(size+i) = 0;
	}
	
// 	if (false)
	{
		while (ifs) {
			ifs.read((char*)&nbytes,8);
			bytes_read+=8;
			if (ifs) {
				if (nbytes!=int8_JCHECK) {
					// read when the check integer was not encountered
	// 				if (myID==0) cout << "reading " << nbytes << " bytes" << endl;
					int64_t * pd = (int64_t*)MyNewCatchChar(__FILE__,__LINE__,nbytes);			// allocate new memory
					if (!pd) {
						cout << "allocation of memory failed!" << endl;
					}
					*(dat+iPointer) = (int64_t) pd;				// store pointer
					*(size+iPointer) = nbytes;					// store size
					ifs.read((char*)(pd),nbytes);				// read data
					if (!ifs) {
						cout << "problem reading from file: " << FileName << endl;
						return (-1);
	// 					cout << "aborting..." << endl;
	// 					exit(-1);
					}
					bytes_read+=nbytes;
					iPointer++;
				}
				else {
	// 				cout << "read " << (double)bytes_read/(double)1000000 << " Mbytes " << endl;
	// 				cout << "check int!" << endl;
					check_int_ok=true;
				}
			}
			if (iPointer==iMax) break; // stop when the desired number of fields was read
		}
	}

	ifs.close();
	if (!check_int_ok) {
		if (iMax>100) {
			cout << "check int missing!" << endl;
		}
	}
	
	return 1;
}

int32_t DS::thread_storage (int32_t n_super_cycle, int32_t myID, int32_t nProcs) {
	int32_t i_store=0;
	int64_t ta=MyGetTime();
	int64_t tb=-1;

	string FileName;
	int64_t n_bytes_loaded=0;
	int64_t n_mol_loaded=0;
	int64_t current_nm_sum=-1;


	for (int32_t i_part=0;i_part<n_part;i_part++) {
		int64_t pdat [16], dat_size [16];

		if (i_part % nProcs == myID) {
			// load part from file
			// cout << "load_"<<myID << "_" << i_part << endl;

			FileName="rho_0.5_T_1.5_species/";
#if defined case_30
			FileName+="case_30/stream_p_";
#elif defined case_63
			FileName+="case_63/stream_p_";
#elif defined case_136
			FileName+="case_136/stream_p_";
#elif defined case_293
			FileName+="case_293/stream_p_";
#elif defined case_630
			FileName+="case_630/stream_p_";
#elif defined case_900
			FileName+="case_900/stream_p_";
#else
			FileName="no_file";
#endif

			FileName+=to_string(i_part);
			cout << FileName << endl;
			int32_t ret=FileFieldsToMem((int64_t*)&pdat,(char*)FileName.c_str(),1,dat_size);

			block_check((double*)pdat[0],1);
			// block_mark((double*)pdat[0]);
			n_mol_loaded+=block_get_nm((double*)pdat[0]);

			// send block to rank 0
			int32_t my_tag=i_part*1000+968;
			int32_t res=MPI_Send((const void*)pdat[0],dat_size[0]/sizeof(int32_t),MPI_INT,0,my_tag,MPI_COMM_WORLD);
			if (dat_size[0]%sizeof(int32_t) !=0 ) {
				cout << "size must be a multiple of sizeof(int32_t)!!!" << endl;
			}

			// cudaMemcpy((void*)((int32_t*)d_in[i_store]),(const void*)pdat[0],dat_size[0],cudaMemcpyHostToDevice);
			delete [] (int64_t*)pdat[0];

			if (ret==-1) {
				cout << "error_reading_file" << endl;
				// the program is in a broken state now... exit (todo)
			}
			else {
				n_bytes_loaded+=dat_size[0];
			}
		}

	}

	return 0;
}
