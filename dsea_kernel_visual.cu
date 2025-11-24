// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include <stdio.h>		// printf
// #include <cub/cub.cuh>
#include <fstream>

using namespace :: std;

__global__ void prepare_visual_rectilinear(double * __restrict__ p_in, double * __restrict__ p_out) {

	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;
	// int32_t n_threads = blockDim.x*gridDim.x;

	int32_t * p_in_i32 = (int32_t *)p_in;
	int64_t * p_in_i64 = (int64_t*)p_in;


	int64_t i_part=p_in_i64[block_i_general_ipart];
	// int32_t n_mol=p_in_i32[block_i_general_nmol*2+0];

	char * p_out_char = (char *)p_out;
	float * p_out_float=(float*)p_out;
	// int * p_out_i32=(int*)p_out;

	// if (global_id==0) {
	// 	p_out_i32[0]=n_mol;
	// 	p_out_i32[1]=i_part;
	// }
	if (global_id<my_n_part*my_n_part) {
		int32_t i_cell=global_id;
		int32_t i_x=i_part;
		int32_t i_y=i_cell/my_n_part;
		int32_t i_z=i_cell-i_y*my_n_part;

		// // load data
		// double px=p_in[block_offset_mol+i_mol*block_doubles_per_mol+0];
		// double py=p_in[block_offset_mol+i_mol*block_doubles_per_mol+1];
		// double pz=p_in[block_offset_mol+i_mol*block_doubles_per_mol+2];

		// // convert to float
		// float f_px=px;
		// float f_py=py;
		// float f_pz=pz;

		// float f_data=0;

		// p_out_float[2+i_mol]=f_data;

		// printf("i_part_%i_%i_%i_%i_\n",i_x,i_y,i_z,i_part);

		int n_mol_cc=p_in_i32[block_offset_nm_cell*2+i_cell];


		int n_mol_species[n_species_max];
		for (int i_species=0;i_species<n_species_max;i_species++) {
			n_mol_species[i_species]=0;
		}

		for (int i_mol=0;i_mol<n_mol_cc;i_mol++) {
			int j_mol=p_in_i32[block_offset_cell_list*2+i_cell*c_mol_max+i_mol];
			// printf("i_%i_%i_%i\n",i_mol,j_mol,n_mol);

			int64_t i_species=p_in_i64[block_offset_mol+j_mol*block_doubles_per_mol+12];
			n_mol_species[i_species]++;
		}

		// if (global_id==0) {
			// printf("%i_%i\n",n_mol_species[0],n_mol_species[1]);
		// }

		// density
		// p_out_float[0*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=n_mol_cc;
		// p_out_char[0*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=n_mol_cc;

		for (int i_species=0;i_species<n_species_max;i_species++) {
			p_out_char[i_species*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=n_mol_species[i_species];
		}


		// double vx_sum=0;
		// double vy_sum=0;
		// double vz_sum=0;
		// double vxvx_sum=0;
		// double vyvy_sum=0;
		// double vzvz_sum=0;
		// for (int i_mol=0;i_mol<n_mol_cc;i_mol++) {
		// 	int j_mol=p_in_i32[block_offset_cell_list*2+i_cell*c_mol_max+i_mol];
		// 	// printf("i_%i_%i_%i\n",i_mol,j_mol,n_mol);

		// 	double vx=p_in[block_offset_mol+j_mol*block_doubles_per_mol+3];
		// 	double vy=p_in[block_offset_mol+j_mol*block_doubles_per_mol+4];
		// 	double vz=p_in[block_offset_mol+j_mol*block_doubles_per_mol+5];

		// 	vx_sum+=vx;
		// 	vy_sum+=vy;
		// 	vz_sum+=vz;
		// 	vxvx_sum+=vx*vx;
		// 	vyvy_sum+=vy*vy;
		// 	vzvz_sum+=vz*vz;
		// }
		// p_out_float[1*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=vx_sum;
		// p_out_float[2*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=vy_sum;
		// p_out_float[3*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=vz_sum;
		// p_out_float[4*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=vxvx_sum;
		// p_out_float[5*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=vyvy_sum;
		// p_out_float[6*my_n_part*my_n_part*my_n_part+i_z*my_n_part*my_n_part+i_y*my_n_part+i_x]=vzvz_sum;
	}

}

// void DS::write_vtk_rectilinear (float * p_data, int32_t n_mol, int32_t i_part, int32_t i_cycle) {
// 	string FileName;
// 	FileName.append("visual/visual_");
// 	FileName+=to_string(i_cycle);
// 	// FileName.append("/visual_");
// 	// FileName+=to_string(i_part);
// 	FileName.append(".vtk");
// 	// cout << "write_vtk_rectilinear" << endl;
// 	ofstream ofs;
// 	ofs.open(FileName, ios::out | ios::binary);
// 	if (ofs) {
// 		ofs << "# vtk DataFile Version 3.0" << endl;
// 		ofs << "vtk output" << endl;
// 		ofs << "ASCII" << endl;
// 		ofs << "DATASET RECTILINEAR_GRID" << endl;
// 		ofs << "DIMENSIONS " << my_n_part+1 << " " << my_n_part+1 << " " << my_n_part+1 << endl;
// 		ofs << "X_COORDINATES " << my_n_part+1 << " float" << endl;
// 		for (int i=0;i<my_n_part+1;i++) {
// 			ofs << i << " ";
// 		}
// 		ofs << endl;
// 		ofs << "Y_COORDINATES " << my_n_part+1 << " float" << endl;
// 		for (int i=0;i<my_n_part+1;i++) {
// 			ofs << i << " ";
// 		}
// 		ofs << endl;
// 		ofs << "Z_COORDINATES " << my_n_part+1 << " float" << endl;
// 		for (int i=0;i<my_n_part+1;i++) {
// 			ofs << i << " ";
// 		}
// 		ofs << endl;

// 		int64_t n_cell_output=my_n_part*my_n_part*my_n_part;

// 		ofs << "CELL_DATA " << n_cell_output << endl;

// 		ofs << "SCALARS density float 1" << endl;
// 		ofs << "LOOKUP_TABLE default" << endl;
// 		for (int64_t i=0;i<n_cell_output;i++) {
// 			ofs << p_data[i] << " ";
// 		}
// 		ofs << endl;

// 		ofs << "SCALARS vx float 1" << endl;
// 		ofs << "LOOKUP_TABLE default" << endl;
// 		for (int64_t i=0;i<n_cell_output;i++) {
// 			ofs << p_data[1*my_n_part*my_n_part*my_n_part+i] << " ";
// 		}
// 		ofs << endl;

// 		ofs << "SCALARS vy float 1" << endl;
// 		ofs << "LOOKUP_TABLE default" << endl;
// 		for (int64_t i=0;i<n_cell_output;i++) {
// 			ofs << p_data[2*my_n_part*my_n_part*my_n_part+i] << " ";
// 		}
// 		ofs << endl;

// 		ofs << "SCALARS vz float 1" << endl;
// 		ofs << "LOOKUP_TABLE default" << endl;
// 		for (int64_t i=0;i<n_cell_output;i++) {
// 			ofs << p_data[3*my_n_part*my_n_part*my_n_part+i] << " ";
// 		}
// 		ofs << endl;

// 		ofs << "SCALARS vxvx float 1" << endl;
// 		ofs << "LOOKUP_TABLE default" << endl;
// 		for (int64_t i=0;i<n_cell_output;i++) {
// 			ofs << p_data[4*my_n_part*my_n_part*my_n_part+i] << " ";
// 		}
// 		ofs << endl;

// 		ofs << "SCALARS vyvy float 1" << endl;
// 		ofs << "LOOKUP_TABLE default" << endl;
// 		for (int64_t i=0;i<n_cell_output;i++) {
// 			ofs << p_data[5*my_n_part*my_n_part*my_n_part+i] << " ";
// 		}
// 		ofs << endl;

// 		ofs << "SCALARS vzvz float 1" << endl;
// 		ofs << "LOOKUP_TABLE default" << endl;
// 		for (int64_t i=0;i<n_cell_output;i++) {
// 			ofs << p_data[6*my_n_part*my_n_part*my_n_part+i] << " ";
// 		}
// 		ofs << endl;
// 		ofs.close();
// 	}
// }

void DS::write_vtr (float * p_data, int32_t i_part, int32_t i_cycle) {
	string FileName;
	FileName.append("visual/visual_");
	FileName+=to_string(i_cycle);
	// FileName.append("/visual_");
	// FileName+=to_string(i_part);
	FileName.append(".vtr");

	ofstream ofs;
	ofs.open(FileName, ios::out | ios::binary);
	if (ofs) {
		int64_t append_offset=0;
		ofs << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		ofs << "<RectilinearGrid WholeExtent=\"" << "0 " << my_n_part << " 0 " << my_n_part << " 0 " << my_n_part << "\">" << endl;
		ofs << "<Piece Extent=\"" << "0 " << my_n_part << " 0 " << my_n_part << " 0 " << my_n_part << "\">" << endl;

		ofs << "<CellData Scalars=\"\" Name=\"a\">";

		// ofs << "<DataArray type=\"Int8\" Name=\"density\" NumberOfComponents=\"1\" format=\"appended\" offset=\"";
		// ofs << append_offset;
		// ofs << "\">";
		// ofs << "</DataArray>" << endl;
		// append_offset+=(my_n_part*block_ncc)*sizeof(char)+sizeof(int64_t);

		for (int i_species=0;i_species<n_species_max;i_species++) {
			ofs << "<DataArray type=\"Int8\" Name=\"density_species_" << i_species << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"";
			ofs << append_offset;
			ofs << "\">";
			ofs << "</DataArray>" << endl;
			append_offset+=(my_n_part*block_ncc)*sizeof(char)+sizeof(int64_t);
		}

		ofs << "</CellData>" << endl;

		ofs << "<Coordinates>" << endl;
		ofs << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"1\" format=\"appended\" offset=\"";
		ofs << append_offset;
		ofs << "\">";
		ofs << "</DataArray>" << endl;
		append_offset+=(my_n_part+1)*sizeof(float)+sizeof(int64_t);

		ofs << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"1\" format=\"appended\" offset=\"";
		ofs << append_offset;
		ofs << "\">";
		// ofs << "\" RangeMin=\"0\" RangeMax=\"1.0\">" << endl;
		ofs << "</DataArray>" << endl;
		append_offset+=(my_n_part+1)*sizeof(float)+sizeof(int64_t);

		ofs << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"1\" format=\"appended\" offset=\"";
		ofs << append_offset;
		ofs << "\">";
		// ofs << "\" RangeMin=\"0\" RangeMax=\"1.0\">" << endl;
		ofs << "</DataArray>" << endl;
		append_offset+=(my_n_part+1)*sizeof(float)+sizeof(int64_t);

		ofs << "</Coordinates>" << endl;

		// ofs << "\" NumberOfCells=\"0\">" << endl;
		// ofs << "<PointData Scalars=\"species\">" << endl;
		// ofs << "<DataArray type=\"Float32\" Name=\"species\" format=\"appended\" offset=\"0\" RangeMin=\"0\" RangeMax=\"6\">" << endl;
		// ofs << "</DataArray>" << endl;
		// ofs << "</PointData>" << endl;
		// ofs << "<Points>" << endl;
		// ofs << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"";
		// ofs << n_mol*sizeof(float)+8;
		// ofs << "\" RangeMin=\"0\" RangeMax=\"1.0\">" << endl;
		// ofs << "</DataArray>" << endl;
		// ofs << "</Points>" << endl;
		// ofs << "<Cells>" << endl;
		// ofs << "<DataArray type=\"Int32\" Name=\"connectivity\"></DataArray>" << endl;
		// ofs << "<DataArray type=\"Int32\" Name=\"offsets\"></DataArray>" << endl;
		// ofs << "<DataArray type=\"UInt8\" Name=\"types\"></DataArray>" << endl;
		// ofs << "</Cells>" << endl;
		ofs << "</Piece>" << endl;
		ofs << "</RectilinearGrid>" << endl;
		ofs << "<AppendedData encoding=\"raw\">" << endl;
		ofs << "_";	// mark start of appended data
		ofs.close();
	}

	// write appended data
	int64_t size_append=0;

	// cell data

	char * p_data_to_write = (char*)p_data;

	for (int i_species=0;i_species<n_species_max;i_species++) {
		size_append=(my_n_part*block_ncc)*sizeof(char);
		MemToFile(&size_append,sizeof(int64_t),(char*)FileName.c_str(),0);
		MemToFile((int64_t*)p_data_to_write,size_append,(char*)FileName.c_str(),0);
		p_data_to_write+=my_n_part*my_n_part*my_n_part;
	}



	// coordinates - same for x,y,z
	float * x_coordinates=new float [my_n_part+1];
	for (int i=0;i<my_n_part+1;i++) {
		x_coordinates[i]=i;
	}
	size_append=(my_n_part+1)*sizeof(float);
	MemToFile(&size_append,sizeof(int64_t),(char*)FileName.c_str(),0);
	MemToFile((int64_t*)x_coordinates,size_append,(char*)FileName.c_str(),0);

	size_append=(my_n_part+1)*sizeof(float);
	MemToFile(&size_append,sizeof(int64_t),(char*)FileName.c_str(),0);
	MemToFile((int64_t*)x_coordinates,size_append,(char*)FileName.c_str(),0);

	size_append=(my_n_part+1)*sizeof(float);
	MemToFile(&size_append,sizeof(int64_t),(char*)FileName.c_str(),0);
	MemToFile((int64_t*)x_coordinates,size_append,(char*)FileName.c_str(),0);
	delete [] x_coordinates;

	// write closing tags
	ofs.open(FileName, ios::out | ios::binary | ios_base::app);
	if (ofs) {
		ofs << "</AppendedData>" << endl;
		ofs << "</VTKFile>" << endl;
		ofs.close();
	}

	// if (i_part==n_part-1) {
	// 	// write pvtu

	// 	FileName.clear();
	// 	FileName.append("visual/visual_");
	// 	FileName+=to_string(i_cycle);
	// 	FileName.append(".pvtu");
	// 	ofs.open(FileName, ios::out | ios::binary);
	// 	if (ofs) {
	// 		ofs << "<?xml version=\"1.0\"?>" << endl;
	// 		ofs << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
	// 		ofs << "<PUnstructuredGrid GhostLevel=\"0\">" << endl;
	// 		ofs << "<PPointData Scalars=\"species\">" << endl;
	// 		ofs << "<PDataArray type=\"Float32\" Name=\"species\"/>" << endl;
	// 		ofs << "</PPointData>" << endl;
	// 		ofs << "<PPoints>" << endl;
	// 		ofs << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << endl;
	// 		ofs << "</PPoints>" << endl;
	// 		for (int32_t i_p=0;i_p<n_part;i_p++) {
	// 			ofs << "<Piece Source=\"";
	// 			ofs << "visual_" << i_cycle << "/visual_" << i_p;
	// 			ofs << ".vtu\"/>" << endl;
	// 		}
	// 		ofs << "</PUnstructuredGrid>" << endl;
	// 		ofs << "</VTKFile>" << endl;

	// 		ofs.close();
	// 	}
	// 	else {
	// 		cout << "error opening file " << FileName << endl;
	// 	}
	// }
}

void DS::caller_output_vtk_rectilinear (double * p_in, double * p_out, cudaStream_t * stream, int32_t threads_per_block, int32_t blockSize, int32_t myID, int32_t i_cycle, int32_t i_part) {
	prepare_visual_rectilinear <<<blockSize,threads_per_block,0,*stream>>> (p_in,p_out);

	// int32_t * p_my_vis_i32=(int32_t*)p_my_vis;
	// float * p_my_vis_float=(float*)p_my_vis;
	if (i_part==(my_n_part-1)) {
		// last part
		float * p_my_vis_float=new float[16*my_n_part*block_ncc];

		cudaDeviceSynchronize();        cudaCheckError(__LINE__,__FILE__);

		size_t copy_size=1;
		copy_size*=my_n_part;
		copy_size*=block_ncc;
		copy_size*=n_species_max;
		copy_size*=sizeof(char);

		// cout << copy_size << endl;
		cudaMemcpy((void*)p_my_vis_float,(const void*)p_out,copy_size,cudaMemcpyDeviceToHost);

	// int32_t n_mol=p_my_vis_i32[0];
	// int32_t i_part=p_my_vis_i32[1];
	// // cout << n_mol << "_" << i_part << endl;

		// string new_dir;
		// new_dir.append("visual/visual_");
		// new_dir+=to_string(i_cycle);

		// boost::filesystem::create_directory(new_dir.c_str());
		write_vtr(p_my_vis_float,0,i_cycle);
		delete [] p_my_vis_float;
	}
}