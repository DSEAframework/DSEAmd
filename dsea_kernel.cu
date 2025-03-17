// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include <stdio.h>		// printf
#include <cub/cub.cuh>
#include <climits>       // for INT_MAX

using namespace :: std;

// print variadic template values
// overload
template<typename T>
void myprint(T head)
{
    std::cout << head << std::endl;
}
// base case: used when pack is non-empty
template<typename T, typename... Ts>
void myprint(T head, Ts... tail)
{
    std::cout << head << std::endl;
    myprint(tail...);
}

// init block, i.e. reset counters, lists...
__global__ void init_block(double * __restrict__ p_in) {
	int32_t * p_in_i32 = (int32_t*) p_in;
	// int64_t * p_in_i64 = (int64_t*) p_in;

	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;

	if (global_id==0) {
		// int64_t i_part=p_in_i64[block_i_general_ipart];
		// printf("init_%ld\n",p_in_i32[block_i_general_nmol*2+0]);

		p_in_i32[block_i_general_nmol*2+0]=0;

		p_in[block_i_general_add_U]=0;
		p_in[block_i_general_add_V]=0;
	}

	if (global_id==0) {
		int32_t n_threads = blockDim.x*gridDim.x;
		if (n_threads<block_ncc) printf("not_enough_threads_X_%i_%i\n",n_threads,block_ncc);
	}

	// molecules per cell
	if (global_id<block_ncc) {
		p_in_i32[block_offset_nm_cell*2+global_id]=0;
	}
}


__global__ void md_v3a(int32_t i_worker, int32_t order_in, int32_t order_out, double * p_sum,
						/*order 0*/ double * __restrict__ p_in_c,
						/*order 1*/ double * __restrict__ p_in_l, double * __restrict__ p_in_r,
						/*order 0*/ double * __restrict__ p_out_c,
						/*order 1*/ double * __restrict__ p_out_l, double * __restrict__ p_out_r,
						int32_t * __restrict__ p_mol_work_list,
						int32_t * __restrict__ p_debug,
						double * __restrict__ p_part_U, double * __restrict__ p_part_V,
						double * __restrict__ p_part_r_U, double * __restrict__ p_part_r_V,
						double * __restrict__ p_part_vel,
						int64_t * __restrict__ p_sum_N,
						double * __restrict__ p_tmp_f) {
	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;
	// int32_t n_threads = blockDim.x*gridDim.x;

	int32_t * p_in_c_i32 = (int32_t *)p_in_c;
	int64_t * p_in_c_i64 = (int64_t*)p_in_c;

	// int32_t * p_in_l_i32 = (int32_t *)p_in_l;

	int32_t * p_in_r_i32 = (int32_t *)p_in_r;

	int64_t i_part=p_in_c_i64[block_i_general_ipart];
	int32_t n_mol=p_in_c_i32[block_i_general_nmol*2+0];

	int32_t grid_nx=p_in_c_i64[block_i_general_nx];
	int32_t grid_ny=p_in_c_i64[block_i_general_ny];
	int32_t grid_nz=p_in_c_i64[block_i_general_nz];

	int32_t local_id = threadIdx.x;

	if (global_id==0) {
		int32_t n_threads = blockDim.x*gridDim.x;
		p_sum_N[i_part]+=n_mol;
		if (n_threads<n_mol) printf("not_enough_threads_%i_%i\n",n_threads,n_mol);
		if (n_threads>2*n_mol) printf("too_many_threads_%i_%i\n",n_threads,n_mol);
	}

	if (global_id<n_mol) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// load molecule
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int32_t list_i_c=p_mol_work_list[global_id*2+0];
		int32_t list_i_mol=p_mol_work_list[global_id*2+1];

		int32_t i_mol=p_in_c_i32[block_offset_cell_list*2+list_i_c*c_mol_max+list_i_mol];

		double px=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+0];
		double py=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+1];
		double pz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+2];
		// double vx=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+3];
		// double vy=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+4];
		// double vz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+5];
		// double fx_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+6];
		// double fy_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+7];
		// double fz_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+8];
		double fx=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+9];
		double fy=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+10];
		double fz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+11];

		int32_t i_c_part=list_i_c;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calcPotForceLJ
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double HILF_3 = md_rc2/md_box2;
		double UCONST=0.004079223;

		double U=0;
		double V=0;
		double Ur=0;
		double Vr=0;

		int ic=i_c_part;

		int ix=i_part;
		int iy=ic/grid_nz;
		int iz=ic-iy*grid_nz;

		// loop over surrounding grid cells
		for (int cx=0;cx<2;cx++) {
		for (int cy=-1;cy<2;cy++) {
		for (int cz=-1;cz<2;cz++) {

			int test_x=ix+cx;
			int test_y=iy+cy;
			int test_z=iz+cz;
			if (test_y<0) test_y+=grid_ny;
			if (test_z<0) test_z+=grid_nz;
			if (test_y>grid_ny-1) test_y-=grid_ny;
			if (test_z>grid_nz-1) test_z-=grid_nz;

			bool do_test=true;
			if (test_x<0) do_test=false;
			if (test_x>grid_nx-1) do_test=false;


			if (do_test==true) {
			// atomicAdd((int*)&p_debug[4],1);

				int ic_test=test_y*grid_nz+test_z;	// cell to test next
				int n_c_test=-1;
				
				if (cx==0) {
					// if ((cy>=0) &&(cz>=0)) {
					// if (cz>=0) {
						n_c_test=p_in_c_i32[block_offset_nm_cell*2+ic_test];		// number of molecules in cell to test
					// }
				}

				// }
				else if (cx==1) {
					n_c_test=p_in_r_i32[block_offset_nm_cell*2+ic_test];		// number of molecules in cell to test
				}

				for (int k=0;k<n_c_test;k++) {		// loop over molecules in test cell
					// int j_mol=p_in_c_i32[block_offset_cell_list*2+ic_test*c_mol_max+k];

					double j_px=0;
					double j_py=0;
					double j_pz=0;
					bool do_comp=true;

					if (cx==0) {
						int j_mol=p_in_c_i32[block_offset_cell_list*2+ic_test*c_mol_max+k];
						j_px=p_in_c[block_offset_mol+j_mol*block_doubles_per_mol+0];
						j_py=p_in_c[block_offset_mol+j_mol*block_doubles_per_mol+1];
						j_pz=p_in_c[block_offset_mol+j_mol*block_doubles_per_mol+2];
						if (i_mol==j_mol) do_comp=false;
						if (j_mol>i_mol) do_comp=false;

						if (do_comp) {
							double DX = px - j_px;
							double DY = py - j_py;
							double DZ = pz - j_pz;

							// DX = DX - round(DX);  // Minimum Image Convention -not needed in x-direction due to specular reflection

							if ((iy==0)||(iy==grid_ny-1)||(iz==0)||(iz==grid_nz-1)) {
							DY = DY - round(DY);  // Minimum Image Convention
							DZ = DZ - round(DZ);  // Minimum Image Convention
							}

							double RIJ2 = DX*DX + DY*DY + DZ*DZ;
							if (RIJ2 <= HILF_3) {
								// atomicAdd((int*)&p_debug[1],1);
								double RIJ2I = 1./(RIJ2*md_box2);
								double RIJ6I = RIJ2I*RIJ2I*RIJ2I;
								double RIJ12I = RIJ6I*RIJ6I;
								double FF = 24.*(2.*RIJ12I-RIJ6I)*RIJ2I*md_box;

								double FIJ_x = FF*DX;
								double FIJ_y = FF*DY;
								double FIJ_z = FF*DZ;

								// actio - current molecule
								fx+=FIJ_x;
								fy+=FIJ_y;
								fz+=FIJ_z;

								// U += (RIJ12I-RIJ6I+UCONST)/2.0;		// half energy
								// V += (2.0*RIJ12I-RIJ6I)/2.0;		// half energy

								//reactio - second molecule in c block

								U += (RIJ12I-RIJ6I+UCONST);		// full energy
								V += (2.0*RIJ12I-RIJ6I);		// full energy

								atomicAdd(&p_tmp_f[j_mol*3+0],-FIJ_x);
								atomicAdd(&p_tmp_f[j_mol*3+1],-FIJ_y);
								atomicAdd(&p_tmp_f[j_mol*3+2],-FIJ_z);

							}
						}

					}
					// else if (cx==-1) {
					// 	int j_mol=p_in_l_i32[block_offset_cell_list*2+ic_test*c_mol_max+k];
					// 	j_px=p_in_l[block_offset_mol+j_mol*block_doubles_per_mol+0];
					// 	j_py=p_in_l[block_offset_mol+j_mol*block_doubles_per_mol+1];
					// 	j_pz=p_in_l[block_offset_mol+j_mol*block_doubles_per_mol+2];
					// }
					else if (cx==1) {
						// 2nd molecule in right block
						int j_mol=p_in_r_i32[block_offset_cell_list*2+ic_test*c_mol_max+k];
						j_px=p_in_r[block_offset_mol+j_mol*block_doubles_per_mol+0];
						j_py=p_in_r[block_offset_mol+j_mol*block_doubles_per_mol+1];
						j_pz=p_in_r[block_offset_mol+j_mol*block_doubles_per_mol+2];

						double DX = px - j_px;
						double DY = py - j_py;
						double DZ = pz - j_pz;

						// DX = DX - round(DX);  // Minimum Image Convention -not needed in x-direction due to specular reflection

						if ((iy==0)||(iy==grid_ny-1)||(iz==0)||(iz==grid_nz-1)) {
							DY = DY - round(DY);  // Minimum Image Convention
							DZ = DZ - round(DZ);  // Minimum Image Convention
						}

						double RIJ2 = DX*DX + DY*DY + DZ*DZ;
						if (RIJ2 <= HILF_3) {
							// atomicAdd((int*)&p_debug[1],1);
							double RIJ2I = 1./(RIJ2*md_box2);
							double RIJ6I = RIJ2I*RIJ2I*RIJ2I;
							double RIJ12I = RIJ6I*RIJ6I;
							double FF = 24.*(2.*RIJ12I-RIJ6I)*RIJ2I*md_box;

							double FIJ_x = FF*DX;
							double FIJ_y = FF*DY;
							double FIJ_z = FF*DZ;

							// actio - current molecule
							fx+=FIJ_x;
							fy+=FIJ_y;
							fz+=FIJ_z;

							U += (RIJ12I-RIJ6I+UCONST)/2.0;		// half energy
							V += (2.0*RIJ12I-RIJ6I)/2.0;		// half energy

							//reactio - second molecule in r block
							atomicAdd(&p_in_r[block_offset_mol+j_mol*block_doubles_per_mol+9],-FIJ_x);
							atomicAdd(&p_in_r[block_offset_mol+j_mol*block_doubles_per_mol+10],-FIJ_y);
							atomicAdd(&p_in_r[block_offset_mol+j_mol*block_doubles_per_mol+11],-FIJ_z);

							Ur += (RIJ12I-RIJ6I+UCONST)/2.0;	// half energy
							Vr += (2.0*RIJ12I-RIJ6I)/2.0;		// half energy

						}
					}
				}
			}
		}
		}
		}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2nd part of verlet integrator: new velocity
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// double HILF_2 = 0.5*md_dt/md_box;
		// vx+=(fx+fx_alt)*HILF_2;
		// vy+=(fy+fy_alt)*HILF_2;
		// vz+=(fz+fz_alt)*HILF_2;

		// double sum_vel=vx*vx+vy*vy+vz*vz;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// store molecule
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+0]=px;
		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+1]=py;
		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+2]=pz;
		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+3]=vx;
		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+4]=vy;
		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+5]=vz;
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+9]=fx;	// f becomes f_alt
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+10]=fy;
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+11]=fz;
	
		p_part_U[global_id]=U;
		p_part_V[global_id]=V;
		// p_part_vel[global_id]=sum_vel;
		p_part_r_U[global_id]=Ur;
		p_part_r_V[global_id]=Vr;
	}
	else {
		// thread does not process a particle
		p_part_U[global_id]=0;
		p_part_V[global_id]=0;
		p_part_vel[global_id]=0;
		p_part_r_U[global_id]=0;
		p_part_r_V[global_id]=0;
	}
}



// __global__ void md_v3aa(int32_t i_worker, int32_t order_in, int32_t order_out, double * p_sum,
// 						/*order 0*/ double * __restrict__ p_in_c,
// 						/*order 1*/ double * __restrict__ p_in_l, double * __restrict__ p_in_r,
// 						/*order 0*/ double * __restrict__ p_out_c,
// 						/*order 1*/ double * __restrict__ p_out_l, double * __restrict__ p_out_r,
// 						int32_t * __restrict__ p_mol_work_list,
// 						int32_t * __restrict__ p_debug,
// 						double * __restrict__ p_part_U, double * __restrict__ p_part_V,
// 						double * __restrict__ p_part_r_U, double * __restrict__ p_part_r_V,
// 						double * __restrict__ p_part_vel,
// 						int64_t * __restrict__ p_sum_N,
// 						double * __restrict__ p_tmp_f) {
// 	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;
// 	// int32_t n_threads = blockDim.x*gridDim.x;

// 	int32_t * p_in_c_i32 = (int32_t *)p_in_c;
// 	// int64_t * p_in_c_i64 = (int64_t*)p_in_c;

// 	// int32_t * p_in_l_i32 = (int32_t *)p_in_l;

// 	// int32_t * p_in_r_i32 = (int32_t *)p_in_r;

// 	// int64_t i_part=p_in_c_i64[block_i_general_ipart];
// 	int32_t n_mol=p_in_c_i32[block_i_general_nmol*2+0];

// 	// int32_t grid_nx=p_in_c_i64[block_i_general_nx];
// 	// int32_t grid_ny=p_in_c_i64[block_i_general_ny];
// 	// int32_t grid_nz=p_in_c_i64[block_i_general_nz];

// 	// int32_t local_id = threadIdx.x;


// 	if (global_id<n_mol) {
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// // load molecule
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 		int32_t list_i_c=p_mol_work_list[global_id*2+0];
// 		int32_t list_i_mol=p_mol_work_list[global_id*2+1];

// 		int32_t i_mol=p_in_c_i32[block_offset_cell_list*2+list_i_c*c_mol_max+list_i_mol];

// 		// double px=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+0];
// 		// double py=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+1];
// 		// double pz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+2];
// 		double vx=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+3];
// 		double vy=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+4];
// 		double vz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+5];
// 		double fx_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+6];
// 		double fy_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+7];
// 		double fz_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+8];
// 		double fx=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+9];
// 		double fy=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+10];
// 		double fz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+11];

// 		double tmp_fx=p_tmp_f[3*i_mol+0];
// 		double tmp_fy=p_tmp_f[3*i_mol+1];
// 		double tmp_fz=p_tmp_f[3*i_mol+2];
// 		fx+=tmp_fx;
// 		fy+=tmp_fy;
// 		fz+=tmp_fz;

// 		p_tmp_f[3*i_mol+0]=0;
// 		p_tmp_f[3*i_mol+1]=0;
// 		p_tmp_f[3*i_mol+2]=0;

// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// // 2nd part of verlet integrator: new velocity
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 		double HILF_2 = 0.5*md_dt/md_box;
// 		vx+=(fx+fx_alt)*HILF_2;
// 		vy+=(fy+fy_alt)*HILF_2;
// 		vz+=(fz+fz_alt)*HILF_2;

// 		double sum_vel=vx*vx+vy*vy+vz*vz;
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// // store molecule
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+0]=px;
// 		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+1]=py;
// 		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+2]=pz;
// 		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+3]=vx;
// 		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+4]=vy;
// 		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+5]=vz;
// 		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+6]=fx;	// f becomes f_alt
// 		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+7]=fy;
// 		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+8]=fz;
	
// 		// p_part_U[global_id]=U;
// 		// p_part_V[global_id]=V;
// 		p_part_vel[global_id]=sum_vel;
// 		// p_part_r_U[global_id]=Ur;
// 		// p_part_r_V[global_id]=Vr;
// 	}
// 	else {
// 		// thread does not process a particle
// 		// p_part_U[global_id]=0;
// 		// p_part_V[global_id]=0;
// 		p_part_vel[global_id]=0;
// 		// p_part_r_U[global_id]=0;
// 		// p_part_r_V[global_id]=0;
// 	}
// }


__global__ void md_v3aa(int32_t i_worker, int32_t order_in, int32_t order_out, double * p_sum,
						/*order 0*/ double * __restrict__ p_in_c,
						/*order 1*/ double * __restrict__ p_in_l, double * __restrict__ p_in_r,
						/*order 0*/ double * __restrict__ p_out_c,
						/*order 1*/ double * __restrict__ p_out_l, double * __restrict__ p_out_r,
						int32_t * __restrict__ p_mol_work_list,
						int32_t * __restrict__ p_debug,
						double * __restrict__ p_part_U, double * __restrict__ p_part_V,
						double * __restrict__ p_part_r_U, double * __restrict__ p_part_r_V,
						double * __restrict__ p_part_vel,
						int64_t * __restrict__ p_sum_N,
						double * __restrict__ p_tmp_f) {
	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;
	// int32_t n_threads = blockDim.x*gridDim.x;

	int32_t * p_in_c_i32 = (int32_t *)p_in_c;
	// int64_t * p_in_c_i64 = (int64_t*)p_in_c;

	// int32_t * p_in_l_i32 = (int32_t *)p_in_l;

	// int32_t * p_in_r_i32 = (int32_t *)p_in_r;

	// int64_t i_part=p_in_c_i64[block_i_general_ipart];
	int32_t n_mol=p_in_c_i32[block_i_general_nmol*2+0];

	// int32_t grid_nx=p_in_c_i64[block_i_general_nx];
	// int32_t grid_ny=p_in_c_i64[block_i_general_ny];
	// int32_t grid_nz=p_in_c_i64[block_i_general_nz];

	// int32_t local_id = threadIdx.x;


	if (global_id<n_mol) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// load molecule
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// int32_t list_i_c=p_mol_work_list[global_id*2+0];
		// int32_t list_i_mol=p_mol_work_list[global_id*2+1];

		int32_t i_mol=global_id;//p_in_c_i32[block_offset_cell_list*2+list_i_c*c_mol_max+list_i_mol];

		// double px=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+0];
		// double py=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+1];
		// double pz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+2];
		double vx=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+3];
		double vy=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+4];
		double vz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+5];
		double fx_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+6];
		double fy_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+7];
		double fz_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+8];
		double fx=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+9];
		double fy=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+10];
		double fz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+11];

		double tmp_fx=p_tmp_f[3*i_mol+0];
		double tmp_fy=p_tmp_f[3*i_mol+1];
		double tmp_fz=p_tmp_f[3*i_mol+2];
		fx+=tmp_fx;
		fy+=tmp_fy;
		fz+=tmp_fz;

		p_tmp_f[3*i_mol+0]=0;
		p_tmp_f[3*i_mol+1]=0;
		p_tmp_f[3*i_mol+2]=0;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2nd part of verlet integrator: new velocity
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double HILF_2 = 0.5*md_dt/md_box;
		vx+=(fx+fx_alt)*HILF_2;
		vy+=(fy+fy_alt)*HILF_2;
		vz+=(fz+fz_alt)*HILF_2;

		double sum_vel=vx*vx+vy*vy+vz*vz;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// store molecule
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+0]=px;
		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+1]=py;
		// p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+2]=pz;
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+3]=vx;
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+4]=vy;
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+5]=vz;
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+6]=fx;	// f becomes f_alt
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+7]=fy;
		p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+8]=fz;
	
		// p_part_U[global_id]=U;
		// p_part_V[global_id]=V;
		p_part_vel[global_id]=sum_vel;
		// p_part_r_U[global_id]=Ur;
		// p_part_r_V[global_id]=Vr;
	}
	else {
		// thread does not process a particle
		// p_part_U[global_id]=0;
		// p_part_V[global_id]=0;
		p_part_vel[global_id]=0;
		// p_part_r_U[global_id]=0;
		// p_part_r_V[global_id]=0;
	}
}

__global__ void md_v3b(int32_t i_worker, int32_t order_in, int32_t order_out, double * p_sum,
						/*order 0*/ double * __restrict__ p_in_c,
						/*order 1*/ double * __restrict__ p_in_l, double * __restrict__ p_in_r,
						/*order 0*/ double * __restrict__ p_out_c,
						/*order 1*/ double * __restrict__ p_out_l, double * __restrict__ p_out_r,
						int32_t * __restrict__ p_mol_work_list,
						int32_t * __restrict__ p_debug						) {
	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;
	// int32_t n_threads = blockDim.x*gridDim.x;

	int32_t * p_in_c_i32 = (int32_t *)p_in_c;
	int64_t * p_in_c_i64 = (int64_t*)p_in_c;

	// int32_t * p_in_l_i32 = (int32_t *)p_in_l;
	// int64_t * p_in_l_i64 = (int64_t*)p_in_l;

	// int32_t * p_in_r_i32 = (int32_t *)p_in_r;
	// int64_t * p_in_r_i64 = (int64_t*)p_in_r;

	int32_t * p_out_c_i32 = (int32_t *)p_out_c;
	// int64_t * p_out_c_i64 = (int64_t*)p_out_c;

	int32_t * p_out_l_i32 = (int32_t *)p_out_l;
	// int64_t * p_out_l_i64 = (int64_t*)p_out_l;

	int32_t * p_out_r_i32 = (int32_t *)p_out_r;
	// int64_t * p_out_r_i64 = (int64_t*)p_out_r;


	double dxdydz=p_in_c[block_i_general_dxdydz];
	// int64_t n_part=p_in_c_i64[block_i_general_npart];
	int64_t i_part=p_in_c_i64[block_i_general_ipart];
	int32_t n_mol=p_in_c_i32[block_i_general_nmol*2+0];

	// int32_t grid_nx=p_in_c_i64[block_i_general_nx];
	// int32_t grid_ny=p_in_c_i64[block_i_general_ny];
	int32_t grid_nz=p_in_c_i64[block_i_general_nz];

	if (global_id==0) {
	// 	int32_t n_threads = blockDim.x*gridDim.x;
		// if (n_threads<n_mol) printf("not_enough_threads_%i_%i\n",n_threads,n_mol);
		// printf("p_sum[0]_%e\n",p_sum[0]);
	}

	// copy meta data to output center
	int32_t general_doubles_to_copy=6;
	if (global_id<general_doubles_to_copy) {
		p_out_c[global_id]=p_in_c[global_id];
	}

	if (global_id<n_mol) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// load molecule
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int32_t list_i_c=p_mol_work_list[global_id*2+0];
		int32_t list_i_mol=p_mol_work_list[global_id*2+1];

		int32_t i_mol=p_in_c_i32[block_offset_cell_list*2+list_i_c*c_mol_max+list_i_mol];

		double px=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+0];
		double py=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+1];
		double pz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+2];
		double vx=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+3];
		double vy=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+4];
		double vz=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+5];
		double fx_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+6];
		double fy_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+7];
		double fz_alt=p_in_c[block_offset_mol+i_mol*block_doubles_per_mol+8];

		// int32_t vdeb=0;

		// int32_t tr=p_in_c_i32[(block_offset_mol+i_mol*block_doubles_per_mol+9)*2+1];
		int32_t i_c_part=list_i_c;

		// thermostat - apply previously determined velocity scale factor
		double TEMP=p_sum[0];

		vx*=TEMP;
		vy*=TEMP;
		vz*=TEMP;

// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// // 1st part of velocity verlet integrator: new position
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		double HILF_1 = 0.5*md_dt*md_dt/md_box;

		double px_new=px+vx*md_dt+fx_alt*HILF_1;
		// specularly reflecting wall in x-direction
		if (px_new<0) {
			// vdeb=1;
			px_new=-px_new;
			vx=-vx;
		}
		else if (px_new>1.0) {
			// vdeb=2;
			double dx=px_new-1.0;
			px_new=1.0-dx;
			vx=-vx;
		}

		double py_new=py+vy*md_dt+fy_alt*HILF_1;
		double pz_new=pz+vz*md_dt+fz_alt*HILF_1;

		// periodic boundaries in y and z-direction
		if (py_new<0) py_new+=1.0;
		else if (py_new>1.0) py_new-=1.0;
		if (pz_new<0) pz_new+=1.0;
		else if (pz_new>1.0) pz_new-=1.0;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// store molecule
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		int32_t icx=px_new/dxdydz;
		int32_t icy=py_new/dxdydz;
		int32_t icz=pz_new/dxdydz;

		i_c_part=icy*grid_nz+icz;

		// if (i_c_part>grid_ny*grid_nz) {
		// 	printf("%i_%i\n",i_c_part,grid_ny*grid_nz);
			// printf("%e_%e_%e\n",px_new,py_new,pz_new);
		// }

		if (icx==i_part) {
			// index of molecule in cell list
			int32_t i_mol_c=atomicAdd((int*)&p_out_c_i32[block_offset_nm_cell*2+i_c_part],1);

			// index of molecule in block list
			int32_t i_store=atomicAdd((int*)&p_out_c_i32[block_i_general_nmol*2+0],1);

			p_out_c_i32[block_offset_cell_list*2+i_c_part*c_mol_max+i_mol_c]=i_store;

			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+0]=px_new;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+1]=py_new;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+2]=pz_new;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+3]=vx;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+4]=vy;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+5]=vz;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+6]=fx_alt;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+7]=fy_alt;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+8]=fz_alt;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+9]=0;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+10]=0;
			p_out_c[block_offset_mol+i_store*block_doubles_per_mol+11]=0;
		}
		else if (icx==i_part+1) {
			// index of molecule in cell list
			int32_t i_mol_c=atomicAdd((int*)&p_out_r_i32[block_offset_nm_cell*2+i_c_part],1);

			// index of molecule in block list
			int32_t i_store=atomicAdd((int*)&p_out_r_i32[block_i_general_nmol*2+0],1);

			p_out_r_i32[block_offset_cell_list*2+i_c_part*c_mol_max+i_mol_c]=i_store;

			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+0]=px_new;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+1]=py_new;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+2]=pz_new;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+3]=vx;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+4]=vy;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+5]=vz;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+6]=fx_alt;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+7]=fy_alt;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+8]=fz_alt;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+9]=0;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+10]=0;
			p_out_r[block_offset_mol+i_store*block_doubles_per_mol+11]=0;
		}
		else if (icx==i_part-1) {
			// index of molecule in cell list
			int32_t i_mol_c=atomicAdd((int*)&p_out_l_i32[block_offset_nm_cell*2+i_c_part],1);

			// index of molecule in block list
			int32_t i_store=atomicAdd((int*)&p_out_l_i32[block_i_general_nmol*2+0],1);

			p_out_l_i32[block_offset_cell_list*2+i_c_part*c_mol_max+i_mol_c]=i_store;

			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+0]=px_new;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+1]=py_new;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+2]=pz_new;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+3]=vx;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+4]=vy;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+5]=vz;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+6]=fx_alt;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+7]=fy_alt;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+8]=fz_alt;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+9]=0;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+10]=0;
			p_out_l[block_offset_mol+i_store*block_doubles_per_mol+11]=0;

		}
		else {
			// molecule crossed neighbouring cells
			// printf("waaaa_%i_%i_%e_%e_%e_%i_%i\n",icx,i_part,px,vx,fx,n_mol,vdeb);
			// printf("waaaa_%i_%i_%e_%e_%e\n",icx,i_part,px_new,py_new,pz_new);
			// printf("waaaa_%i_%i_%e_%e_%e\n",icx,i_part,vx,vy,vz);
		}
	}
}






__global__ void md_thermo_a(double * __restrict__ p_in, double * __restrict__ p_sum) {

	__shared__ double sdata[32];
	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;
	int32_t tid=threadIdx.x;
	int32_t n_threads = blockDim.x*gridDim.x;

	int32_t * p_in_i32 = (int32_t *)p_in;
	// int64_t * p_in_i64 = (int64_t*)p_in;

	int32_t n_mol=p_in_i32[block_i_general_nmol*2+0];
	int32_t n_loop=n_mol/n_threads+1;

	double sum=0;

	if (global_id==0) {
		// printf("n_loop_%i\n",n_loop);
	}

	for (int i_loop=0;i_loop<n_loop;i_loop++) {
		int32_t i_mol=i_loop*n_threads+global_id;
		if (i_mol<n_mol) {
			double vx=p_in[block_offset_mol+i_mol*block_doubles_per_mol+3];
			double vy=p_in[block_offset_mol+i_mol*block_doubles_per_mol+4];
			double vz=p_in[block_offset_mol+i_mol*block_doubles_per_mol+5];
			sum+=vx*vx+vy*vy+vz*vz;
		}
	}

	sdata[tid]=sum;
	__syncthreads();

	if (tid==0) {
		sum=0;
		for (int32_t i=0;i<blockDim.x;i++) {
			sum+=sdata[i];
		}
		p_sum[blockIdx.x]=sum;
	}
}


__global__ void md_thermo_b(double * __restrict__ p_in, double * __restrict__ p_sum) {
	int32_t * p_in_i32 = (int32_t *)p_in;
	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;
	int32_t n_mol=p_in_i32[block_i_general_nmol*2+0];

	double sum=0;
	if (global_id==0) {
		for (int32_t i_block=0;i_block<128;i_block++) {
			sum+=p_sum[i_block];
		}

		p_sum[1]=sum;

		sum *= md_box2;
		double TEMP = sqrt(3.*n_mol*md_T/sum);

		// printf("%e_%e_#TEMP\n",sum,TEMP);
		if (TEMP<1e-5) printf("bad_temp_%e",TEMP);

		p_sum[0]=TEMP;

	}

}

// __global__ void md_debug_a(int32_t * __restrict__ a, int32_t * __restrict__ b, int32_t n) {
// 	for (int i=0;i<n;i++) {
// 		printf("deba_%i_%i_%i\n",i,a[i],b[i]);
// 	}
// }

__global__ void stat_collect(	double * __restrict__ p_sum_u, double * __restrict__ p_sum_v, double * __restrict__ p_sum_vel,
								double * __restrict__ p_result, int32_t i_part, int32_t n_part) {
	// printf("debb_%e_%e_%e\n",a[0],b[0],b[1]);
	p_sum_u[i_part]+=p_result[0];
	p_sum_v[i_part]+=p_result[1];
	p_sum_vel[i_part]+=p_result[2];
	if (i_part<n_part-1) {
		p_sum_u[i_part+1]+=p_result[3];
		p_sum_v[i_part+1]+=p_result[4];
	}
}


// each thread fills the list of molecules for one cell
__global__ void md_fill_mol_work_list(int32_t * __restrict__ mol_list, int32_t * __restrict__ nm, int32_t * __restrict__ prefix_sum, int32_t n) {
	int32_t global_id = blockIdx.x*blockDim.x+threadIdx.x;
	// int32_t block_id=blockIdx.x;
	// int32_t tid=threadIdx.x;
	// int32_t n_mol=p_in_i32[block_i_general_nmol*2+0];

	if (global_id<n) {
		int32_t n_mol=nm[global_id];
		int32_t offset=prefix_sum[global_id];

		for (int i_mol=0;i_mol<n_mol;i_mol++) {
			mol_list[offset*2+i_mol*2+0]=global_id;	// cell
			mol_list[offset*2+i_mol*2+1]=i_mol;		// molecule in cell
		}
	}

	
}

void DS::caller_worker (double ** p_in, double ** p_out, int32_t i_part, int32_t i_super_cycle,
						int32_t order_in, int32_t order_out, int32_t iworker, int32_t nworker,
						cudaStream_t * stream, int32_t threads_per_block, int32_t blockSize, int32_t myID) {

// the order of arrays in p_in and p_out is:
// center, left, right, left-left, right-right, left-left-left, right-right-right, and so on
// entries can be NULL when invalid

// launch kernel:
// MatAdd<<<numBlocks, threadsPerBlock>>>(

	// init sum arrays
	if ((i_part==0)&&(i_super_cycle==0)) {
		cudaDeviceSynchronize();	cudaCheckError(__LINE__,__FILE__);
		cudaMemset((void*)d_sum_u,0,my_n_part*sizeof(double));
		cudaMemset((void*)d_sum_v,0,my_n_part*sizeof(double));
		cudaMemset((void*)d_sum_vel,0,my_n_part*sizeof(double));
		cudaMemset((void*)d_sum_N,0,my_n_part*sizeof(int64_t));
		n_cycle_sample=0;
	}

	if ((i_part==0)&&(i_super_cycle==1000)) {
		// reset sample
		cudaDeviceSynchronize();	cudaCheckError(__LINE__,__FILE__);
		cudaMemset((void*)d_sum_u,0,my_n_part*sizeof(double));
		cudaMemset((void*)d_sum_v,0,my_n_part*sizeof(double));
		cudaMemset((void*)d_sum_vel,0,my_n_part*sizeof(double));
		cudaMemset((void*)d_sum_N,0,my_n_part*sizeof(int64_t));
		n_cycle_sample=0;
	}

	// if (false)
	if ((i_part==0)&&(i_super_cycle%100==0)&&(i_super_cycle>0)) {
	// if (i_part==0) {
		cudaDeviceSynchronize();	cudaCheckError(__LINE__,__FILE__);
		cudaMemcpy((void*)h_sum_u,(const void*)d_sum_u,my_n_part*sizeof(double),cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);
		cudaMemcpy((void*)h_sum_v,(const void*)d_sum_v,my_n_part*sizeof(double),cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);
		cudaMemcpy((void*)h_sum_vel,(const void*)d_sum_vel,my_n_part*sizeof(double),cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);
		cudaMemcpy((void*)h_sum_N,(const void*)d_sum_N,my_n_part*sizeof(int64_t),cudaMemcpyDeviceToHost);		cudaCheckError(__LINE__,__FILE__);

		for (int32_t k=0;k<3*my_n_part;k++) {
			// vdata_send_a[k]=0;
			vdata_rec_a[k]=0;
		}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

		for (int32_t i_part=0;i_part<my_n_part;i_part++) {
			vdata_send_a[i_part*3+0]=h_sum_u[i_part];
			vdata_send_a[i_part*3+1]=h_sum_v[i_part];
			vdata_send_a[i_part*3+2]=h_sum_vel[i_part];
		}

		for (int32_t k=0;k<my_n_part;k++) {
			vdata_send_b[k]=0;
			vdata_rec_b[k]=0;
		}
		for (int32_t i_part=0;i_part<my_n_part;i_part++) {
			vdata_send_b[i_part]=h_sum_N[i_part];
		}

		int32_t mpi_ret=-1;
		mpi_ret=MPI_Ireduce(vdata_send_a,vdata_rec_a,3*my_n_part,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,&request_results[0]);
		if (mpi_ret!=MPI_SUCCESS) cout << "waa" << endl;
		mpi_ret=MPI_Ireduce(vdata_send_b,vdata_rec_b,my_n_part,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD,&request_results[1]);
		if (mpi_ret!=MPI_SUCCESS) cout << "waa" << endl;
		outstanding_results=true;
		outstanding_results_n_cycle_sample=n_cycle_sample*8;
	}


	if (outstanding_results) {
		MPI_Status stat[2];
        int32_t flag=false;
		MPI_Testall(2,request_results,&flag,stat);
		if ((flag==true)&&(myID==0)) {
			outstanding_results=false;

			// output global state
			{
				double sum_u=0;
				double sum_v=0;
				double sum_vel=0;

				int64_t stat_N=0;

				for (int32_t i_part=0;i_part<my_n_part;i_part++) {
					sum_u+=vdata_rec_a[i_part*3+0];
					sum_v+=vdata_rec_a[i_part*3+1];
					sum_vel+=vdata_rec_a[i_part*3+2];
					stat_N+=vdata_rec_b[i_part];
				}

				if (stat_N>0) {
					sum_u*=4.0;
					sum_v*=8.0;

					double SUMV = sum_vel*md_box2;
					double SUMK = 0.5*SUMV;
					double SUMU = sum_u/(double)stat_N;
					double SUMP = sum_v/(double)stat_N;

					double TEMP = 2.0*SUMK/3.0;
					TEMP/=(double)stat_N;
					double EPOT = SUMU;
					double EKIN = SUMK/((double)stat_N);
					double ETOT = EPOT+EKIN;
					double DRUCK = TEMP*md_rho+md_rho*SUMP;

					cout << "#stat_sum: " << i_super_cycle << " " << (double)stat_N << " " << TEMP << " " << md_rho << " " << DRUCK << " " << EPOT << " " << EKIN << " " << ETOT << endl;
				}
			}
			// output state per part

			{
				for (int32_t i_part=0;i_part<my_n_part;i_part++) {
					double sum_u=vdata_rec_a[i_part*3+0];
					double sum_v=vdata_rec_a[i_part*3+1];
					double sum_vel=vdata_rec_a[i_part*3+2];
					int64_t stat_N=vdata_rec_b[i_part];

					if (stat_N!=0) {
						sum_u*=4.0;
						sum_v*=8.0;

						stat_N/=outstanding_results_n_cycle_sample;

						double SUMV = sum_vel*md_box2;
						double SUMK = 0.5*SUMV;
						double SUMU = sum_u/(double)stat_N;
						double SUMP = sum_v/(double)stat_N;

						double TEMP = 2.0*SUMK/3.0;
						TEMP/=(double)stat_N;
						double NUM=outstanding_results_n_cycle_sample;
						TEMP/=NUM;
						double EPOT = SUMU/NUM;
						double EKIN = SUMK/((double)stat_N*NUM);
						double ETOT = EPOT+EKIN;
						double DRUCK = TEMP*md_rho+md_rho*SUMP/NUM;
						double part_rho=(double)stat_N/(md_box2*md_box/((double)my_n_part));

						cout << "#stat_part: " << i_super_cycle << " " << i_part << " " << TEMP << " " << part_rho << " " << DRUCK << " " << EPOT << " " << EKIN << " " << ETOT << endl;
					}
					else {
						cout << "#stat_part: no data" << endl;
					}
				}
			}
		}
	}

	// init output arrays
	if (i_part==0) {
		for (int32_t i=0;i<order_out+1;i++) init_block<<<block_ncc/threads_per_block+1,threads_per_block,0,*stream>>>(p_out[i*2]);
	}
	else if (i_part>=(n_part-order_out)) {
		// no init required

		if (p_out[order_out*2]!=(double*)-1) init_block<<<block_ncc/threads_per_block+1,threads_per_block,0,*stream>>>(p_out[order_out*2]);
	}
	else {
		init_block<<<block_ncc/threads_per_block+1,threads_per_block,0,*stream>>>(p_out[order_out*2]);
	}

	// if (false)
	{	// molecular dynamics
	
	// compute prefix sum for number of molecules
	int32_t * my_d_in=(int32_t*)p_in[0];
	if (size_temp==0) {
		size_t my_size_temp=0;
		// determine size_temp and allocate d_temp
		cub::DeviceScan::ExclusiveSum((void*)NULL,my_size_temp,&my_d_in[block_offset_nm_cell*2],(int32_t*)d_prefix_sum,block_ncc,*stream);
		if (my_size_temp>size_temp) size_temp=my_size_temp;

		my_size_temp=0;
		cub::DeviceReduce::Sum((void*)NULL,my_size_temp,(double*)d_part_u,(double*)d_sum_u,blockSize*threads_per_block,*stream);
		if (my_size_temp>size_temp) size_temp=my_size_temp;

		cudaMalloc((void**)&d_temp,size_temp);			cudaCheckError(__LINE__,__FILE__);


// // L2 prefetching
// size_t dev_res;
// cudaDeviceGetLimit(&dev_res,cudaLimitMaxL2FetchGranularity); cudaCheckError(__LINE__,__FILE__);
// cout << "cudaLimitMaxL2FetchGranularity:" << dev_res << endl;
// dev_res=16;
// cudaDeviceSetLimit(cudaLimitMaxL2FetchGranularity,dev_res); cudaCheckError(__LINE__,__FILE__);
// cudaDeviceGetLimit(&dev_res,cudaLimitMaxL2FetchGranularity); cudaCheckError(__LINE__,__FILE__);
// cout << "cudaLimitMaxL2FetchGranularity:" << dev_res << endl;

	}
	cub::DeviceScan::ExclusiveSum((void*)d_temp,size_temp,&my_d_in[block_offset_nm_cell*2],(int32_t*)d_prefix_sum,block_ncc,*stream);
	// cout << "size_temp_b" << size_temp << endl;

	// fill mol work list
	md_fill_mol_work_list <<<block_ncc/threads_per_block+1,threads_per_block,0,*stream>>>((int32_t*)d_mol_list,&my_d_in[block_offset_nm_cell*2],(int32_t*)d_prefix_sum,block_ncc);


	// if (i_part==0) {
	// 	md_debug_a <<<1,1,0,*stream>>>(&my_d_in[block_offset_nm_cell*2],(int*)d_prefix_sum,block_ncc);
	// 	md_debug_b <<<1,1,0,*stream>>>((int*)d_mol_list,block_ncc*c_mol_max);
	// }

	// V3 does
	// md_v3a
	// cell_list_clear(ncells);
	// cell_list_compute_c(ncells);
	// cell_list_build(ncells);
	// calcPotForceLJ();
	// verlet2();
	// md_v3b
	// getStateValues(NUM);
	// if (FLGENSEM_NVT) { scaleVelocity(); }
	// verlet1();

	md_v3a <<<blockSize,threads_per_block,0,*stream>>>(iworker,order_in,order_out,(double*)d_sum,
												p_in[0],p_in[1],p_in[2],
												p_out[0],p_out[1],p_out[2],
												(int32_t*)d_mol_list,
												(int32_t*)d_debug,
												(double*)d_part_u,(double*)d_part_v,
												(double*)d_part_r_u,(double*)d_part_r_v,
												(double*)d_part_vel,
												(int64_t*)d_sum_N,
												(double*)d_tmp_f);

	md_v3aa <<<blockSize,threads_per_block,0,*stream>>>(iworker,order_in,order_out,(double*)d_sum,
												p_in[0],p_in[1],p_in[2],
												p_out[0],p_out[1],p_out[2],
												(int32_t*)d_mol_list,
												(int32_t*)d_debug,
												(double*)d_part_u,(double*)d_part_v,
												(double*)d_part_r_u,(double*)d_part_r_v,
												(double*)d_part_vel,
												(int64_t*)d_sum_N,
												(double*)d_tmp_f);


	cub::DeviceReduce::Sum((void*)d_temp,size_temp,(double*)d_part_u,(double*)d_result+0,blockSize*threads_per_block,*stream);
	cub::DeviceReduce::Sum((void*)d_temp,size_temp,(double*)d_part_v,(double*)d_result+1,blockSize*threads_per_block,*stream);
	cub::DeviceReduce::Sum((void*)d_temp,size_temp,(double*)d_part_vel,(double*)d_result+2,blockSize*threads_per_block,*stream);
	cub::DeviceReduce::Sum((void*)d_temp,size_temp,(double*)d_part_r_u,(double*)d_result+3,blockSize*threads_per_block,*stream);
	cub::DeviceReduce::Sum((void*)d_temp,size_temp,(double*)d_part_r_v,(double*)d_result+4,blockSize*threads_per_block,*stream);

	// thermostat required for all versions
	md_thermo_a <<<128,threads_per_block,0,*stream>>>(p_in[0],(double*)d_sum);
	md_thermo_b <<<1,threads_per_block,0,*stream>>>(p_in[0],(double*)d_sum);

	stat_collect <<<1,1,0,*stream>>> (	(double*)d_sum_u,(double*)d_sum_v,(double*)d_sum_vel,
										(double*)d_result,i_part,my_n_part);

	md_v3b <<<blockSize,threads_per_block,0,*stream>>>(iworker,order_in,order_out,(double*)d_sum,
												p_in[0],p_in[1],p_in[2],
												p_out[0],p_out[1],p_out[2],
												(int32_t*)d_mol_list,
												(int32_t*)d_debug);
	}

	if (i_part==0) {
		n_cycle_sample++;
	}
}
