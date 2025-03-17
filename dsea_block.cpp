// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

// This routine determines the actual size of a block in bytes. The maximum size is that of the storage, i.e. store_size bytes
int64_t DS::get_block_size_host(double * p_data) {
	int32_t * p_data_i32=(int32_t*)p_data;

	int64_t block_size=block_offset_mol*sizeof(double);

	int32_t block_n_mol=p_data_i32[block_i_general_nmol*2+0];
	// cout << "block_n_mol_" << block_n_mol << endl;
	block_size+=block_doubles_per_mol*block_n_mol*sizeof(double);

	return block_size;
}

// This routine determines the block size of a block stored in GPU memory
int64_t DS::get_block_size_device(double * p_data) {
	int64_t i_tmp=-1;
	cudaMemcpy((void*)&i_tmp,(const void*)&p_data[block_i_general_nmol],sizeof(int64_t),cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);

	int64_t block_size=block_offset_mol*sizeof(double);
	int32_t block_n_mol=i_tmp;
	// cout << "block_n_mol_" << block_n_mol << endl;
	block_size+=block_doubles_per_mol*block_n_mol*sizeof(double);

	return block_size;
}

// This routine determines the block size of a block stored in GPU memory
int64_t DS::get_block_nm_device(double * p_data) {
	int64_t i_tmp=-1;
	cudaMemcpy((void*)&i_tmp,(const void*)&p_data[block_i_general_nmol],sizeof(int64_t),cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);
	return i_tmp;
}

// This routine checks the block for consistency
int32_t DS::block_check(double * p_data, int32_t pos) {
	return 0;
	bool block_error=false;

	int32_t * p_data_i32=(int32_t*)p_data;
	int64_t * p_data_i64=(int64_t*)p_data;
	int32_t * p_nm_cell = (int32_t*) &p_data[block_offset_nm_cell];

	int64_t i_part=p_data_i64[block_i_general_ipart];
	cout << "check_" << i_part << "_@" << pos;

	// check number of molecules is correct
	int32_t nm=p_data_i32[block_i_general_nmol*2+0];
	int32_t nm_sum=0;
	int32_t nm_max=-1;
cout << "nm:" << nm << endl;

	for (int32_t i_c=0;i_c<block_ncc;i_c++) {
		nm_sum+=p_nm_cell[i_c];
		if (p_nm_cell[i_c]>nm_max) nm_max=p_nm_cell[i_c];
	}
	if (nm!=nm_sum) {
		cout << "nm_sum incorrect_" << nm_sum << "_" << nm << endl;
		block_error=true;
		for (int32_t i_c=0;i_c<block_ncc;i_c++) {
			cout << i_c << "_" << p_nm_cell[i_c] << endl;;
		}
	}

	cout << "_" << nm << "_" << nm_max << "_";

	// check for cell list overflow
	for (int32_t i_c=0;i_c<block_ncc;i_c++) {
		if (p_nm_cell[i_c]>c_mol_max) {
			cout << "cell_list_overflow!" << endl;;
			block_error=true;
		}
	}

	// check if molecules are in correct cell
	double v_max=0;

	double dxdydz=p_data[block_i_general_dxdydz];

	if (dxdydz==0) return 0;

	int64_t grid_nz=p_data_i64[block_i_general_nz];
	for (int32_t i_c=0;i_c<block_ncc;i_c++) {
		int32_t * cell_list=&p_data_i32[block_offset_cell_list*2+i_c*c_mol_max];
		for (int32_t i_mol=0;i_mol<p_nm_cell[i_c];i_mol++) {
			// load molecule
			int j_mol=cell_list[i_mol];
			double px=p_data[block_offset_mol+j_mol*block_doubles_per_mol+0];
			double py=p_data[block_offset_mol+j_mol*block_doubles_per_mol+1];
			double pz=p_data[block_offset_mol+j_mol*block_doubles_per_mol+2];
			double vx=p_data[block_offset_mol+j_mol*block_doubles_per_mol+3];
			double vy=p_data[block_offset_mol+j_mol*block_doubles_per_mol+4];
			double vz=p_data[block_offset_mol+j_mol*block_doubles_per_mol+5];
			// double fx=p_data[block_offset_mol+j_mol*block_doubles_per_mol+6];
			// double fy=p_data[block_offset_mol+j_mol*block_doubles_per_mol+7];
			// double fz=p_data[block_offset_mol+j_mol*block_doubles_per_mol+8];
			// double fx_alt=p_data[block_offset_mol+j_mol*block_doubles_per_mol+9];
			// double fy_alt=p_data[block_offset_mol+j_mol*block_doubles_per_mol+10];
			// double fz_alt=p_data[block_offset_mol+j_mol*block_doubles_per_mol+11];

			// int32_t i_c_part=p_data_i32[(block_offset_mol+j_mol*block_doubles_per_mol+9)*2];

			int32_t icx=px/dxdydz;
			int32_t icy=py/dxdydz;
			int32_t icz=pz/dxdydz;

			if (icx!=i_part) {
				cout << "molecule_in_wrong_part_" << icx << "_" << i_part << endl;
				cout << px << " " << dxdydz << endl;
				cout << i_c << " " << i_mol << " " << j_mol << endl;
				block_error=true;
			}

			int32_t i_c_part_check=icy*grid_nz+icz;
			// if (i_c_part_check!=i_c_part) {
			// 	cout << "block_check_wrong_cell_mol_" << i_c_part_check << "_" << i_c_part << endl;
			// 	block_error=true;
			// }
			if (i_c_part_check!=i_c) {
				cout << "block_check_wrong_cell_part" << i_c_part_check << "_" << i_c << " " << icy << "_" << icz << endl;
				cout << dxdydz << endl;
				cout << i_c << "_" << i_mol << "_" << j_mol << endl;
				cout << px << "_" << py << "_" << pz << endl;
				block_error=true;
			}

			bool bounds_prob=false;
			if (px<0) bounds_prob=true;
			if (py<0) bounds_prob=true;
			if (pz<0) bounds_prob=true;
			if (px>1.0) bounds_prob=true;
			if (py>1.0) bounds_prob=true;
			if (pz>1.0) bounds_prob=true;
			if (bounds_prob) {
				cout << "out_of_bounds_"<< px << "_" << py << "_" << pz << endl;
				block_error=true;
			}

			double v=sqrt(vx*vx+vy*vy+vz*vz);
			if (v>v_max) v_max=v;

		}
	}
	cout << "v_max_" << v_max << endl;
	if (block_error==true) cout << "block_contains_errors" << i_part << "_@" << pos << endl;
	return 0;
}

int32_t DS::block_check_device (double* data, int32_t pos) {
	if (data==(double*)-1) return 0;
	cout << "block_check_device disabled" << endl;
	return 0;
	char * my_block = MyNewCatchChar(__FILE__,__LINE__,store_size);
	cudaDeviceSynchronize();	cudaCheckError(__LINE__,__FILE__);
	cudaMemcpy((void*)my_block,(const void*)data,store_size,cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);
	block_check((double*)my_block,pos);
	delete [] my_block;
	return 0;
}


// This routine checks the block for consistency
int32_t DS::block_get_nm(double * p_data) {
	int32_t * p_data_i32=(int32_t*)p_data;
	int64_t * p_data_i64=(int64_t*)p_data;
	int32_t * p_nm_cell = (int32_t*) &p_data[block_offset_nm_cell];

	// check number of molecules is correct
	int32_t nm=p_data_i32[block_i_general_nmol*2+0];

	return nm;
}

int32_t DS::block_mark(double * p_data) {
	int32_t * p_data_i32=(int32_t*)p_data;
	int64_t * p_data_i64=(int64_t*)p_data;
	int32_t * p_nm_cell = (int32_t*) &p_data[block_offset_nm_cell];

	int64_t i_part=p_data_i64[block_i_general_ipart];
	cout << "check_" << i_part;


	for (int32_t i_c=0;i_c<block_ncc;i_c++) {
		int32_t * cell_list=&p_data_i32[block_offset_cell_list*2+i_c*c_mol_max];
		for (int32_t i_mol=0;i_mol<p_nm_cell[i_c];i_mol++) {
			// load molecule
			int j_mol=cell_list[i_mol];
			// double px=p_data[block_offset_mol+j_mol*block_doubles_per_mol+0];
			// double py=p_data[block_offset_mol+j_mol*block_doubles_per_mol+1];
			// double pz=p_data[block_offset_mol+j_mol*block_doubles_per_mol+2];
			// double vx=p_data[block_offset_mol+j_mol*block_doubles_per_mol+3];
			// double vy=p_data[block_offset_mol+j_mol*block_doubles_per_mol+4];
			// double vz=p_data[block_offset_mol+j_mol*block_doubles_per_mol+5];

			// int32_t i_c_part=p_data_i32[(block_offset_mol+j_mol*block_doubles_per_mol+9)*2];
			p_data_i32[(block_offset_mol+j_mol*block_doubles_per_mol+9)*2+1]=0;

		}
	}
	return 0;
}

int32_t DS::block_markcheck(double * p_data) {
	int32_t * p_data_i32=(int32_t*)p_data;
	int64_t * p_data_i64=(int64_t*)p_data;
	int32_t * p_nm_cell = (int32_t*) &p_data[block_offset_nm_cell];

	int64_t i_part=p_data_i64[block_i_general_ipart];
	cout << "check_" << i_part;


	for (int32_t i_c=0;i_c<block_ncc;i_c++) {
		int32_t * cell_list=&p_data_i32[block_offset_cell_list*2+i_c*c_mol_max];
		for (int32_t i_mol=0;i_mol<p_nm_cell[i_c];i_mol++) {
			// load molecule
			int j_mol=cell_list[i_mol];
			// double px=p_data[block_offset_mol+j_mol*block_doubles_per_mol+0];
			// double py=p_data[block_offset_mol+j_mol*block_doubles_per_mol+1];
			// double pz=p_data[block_offset_mol+j_mol*block_doubles_per_mol+2];
			// double vx=p_data[block_offset_mol+j_mol*block_doubles_per_mol+3];
			// double vy=p_data[block_offset_mol+j_mol*block_doubles_per_mol+4];
			// double vz=p_data[block_offset_mol+j_mol*block_doubles_per_mol+5];

			// int32_t i_c_part=p_data_i32[(block_offset_mol+j_mol*block_doubles_per_mol+9)*2];
			int mark=p_data_i32[(block_offset_mol+j_mol*block_doubles_per_mol+9)*2+1];

			if (mark!=1) cout << "waa_" << mark << endl;

		}
	}
	return 0;
}