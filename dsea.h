#include <cstdint>
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#include <mpi.h>

// parameters for algorithm
// #define case_30			// 1e5
//#define case_63			// 1e6
// #define case_136		// 1e7
//#define case_293		// 1e8
// #define case_630		// 1e9
//#define case_900		// 2.9e9

#include "dsea_case.h"

// select test md case
// rho=0.75, T=0.8
/*
#define md_T							0.8			// reduced temperature
#define md_rho							0.75		// reduced density

#if defined case_30
	#define md_box						52.4148
	#define block_ncc					(20*20)
	#define md_N						(30*30*30*4)
 	#define my_n_part					20
#elif defined case_63
	#define md_box						1.1007114e2
	#define block_ncc					(44*44)
	#define md_N						(63*63*63*4)
 	#define my_n_part					44
#elif defined case_136
	#define md_box						237.614
	#define block_ncc					(95*95+1)
	#define md_N						(136*136*136*4)
 	#define my_n_part					95
#elif defined case_293
	#define md_box						511.918
	#define block_ncc					(204*204)
	#define md_N						(293*293*293*4)
 	#define my_n_part					204
#elif defined case_630
	#define md_box						1100.71
	#define block_ncc					(440*440)
	#define md_N						(630*630*630*4)
#else
	#define md_box						-1
	#define block_ncc					(0)
#endif
*/

// rho=0.5, T=1.5
#define md_T							1.5		// reduced temperature
#define md_rho							0.5		// reduced density
#if defined case_30
	#define md_box						60.0
	#define block_ncc					(23*23+1)
	#define md_N						(30*30*30*4)
 	#define my_n_part					23
#elif defined case_63
	#define md_box						126.0
	#define block_ncc					(50*50)
	#define md_N						(63*63*63*4)
 	#define my_n_part					50
#elif defined case_136
	#define md_box						272.0
	#define block_ncc					(108*108)
	#define md_N						(136*136*136*4)
 	#define my_n_part					108
#elif defined case_293
	#define md_box						586.0
	#define block_ncc					(234*234)
	#define md_N						(293*293*293*4)
 	#define my_n_part					234
#elif defined case_630
	#define md_box						1260
	#define block_ncc					(504*504)
	#define md_N						(630*630*630*4)
 	#define my_n_part					504
#elif defined case_900
	#define md_box						1800
	#define block_ncc					(719*719+1)
	#define md_N						(900*900*900*4)
 	#define my_n_part					719
#else
	// #define md_box						-1
	// #define block_ncc					(0)
#endif




#define md_box2							(md_box*md_box)
#define md_rc							2.5
#define md_rc2							(md_rc*md_rc)

#define md_dt							3.0e-3		// time step

// constants block_* must be adapted to the data set
#define block_general_doubles		16
#define block_doubles_per_mol		16

#define c_mol_max					64
#define block_cell_list_ints		(block_ncc*c_mol_max)	// [int32_t]

#define block_i_general_dxdydz		0
#define block_i_general_nx			1
#define block_i_general_ny			2
#define block_i_general_nz			3
#define block_i_general_npart		4
#define block_i_general_ipart		5
#define block_i_general_nmol		6
#define block_i_general_nstep		7
#define block_i_general_add_U		8
#define block_i_general_add_V		9

#define block_offset_nm_cell		block_general_doubles		// [double]
#define block_offset_cell_list		(block_offset_nm_cell+(block_ncc*sizeof(int32_t))/sizeof(double))				// [double]
#define block_offset_mol			(block_offset_cell_list +(block_cell_list_ints*sizeof(int32_t))/sizeof(double))	// [double]


// no modifications required below

#define CHECK_INT					123456

#define mem_state_free			-1
#define mem_state_ready			-2
#define mem_state_ready_b		-3
#define mem_state_bussy			-4

#define mem_type_input			1
#define mem_type_output			2
#define mem_type_store			3
#define mem_type_worker			4

#define n_worker_event			64

#define n_h_store_max			4096

struct mem_info {
	int32_t type;
	int32_t i_part;
	int32_t state;
	int32_t n_use;
	int32_t i_event;
	int32_t i_cycle;
};

class DS {
	public:
		// functions
		// constructor, destructor
		DS(int32_t igpu, int32_t nworker, int32_t npart, int32_t o_in, int32_t o_out, int32_t myid, int32_t nprocs, int32_t nrails);
		~DS();
		int64_t MyGetTime();

		int32_t thread_input(int32_t n_super_cycle, int32_t myID, int32_t nProcs);
		int32_t thread_output(int32_t n_super_cycle, int32_t myID, int32_t nProcs);
		int32_t thread_input_ucx(int argc, char ** argv, int32_t n_super_cycle, int32_t myID, int32_t nProcs);
		int32_t thread_output_ucx(int argc, char ** argv, int32_t n_super_cycle, int32_t myID, int32_t nProcs);
		int32_t thread_main(int32_t n_super_cycle, int32_t order_in, int32_t order_out, int32_t myID);
		int32_t thread_storage (int32_t n_super_cycle, int32_t myID, int32_t nProcs);

		int32_t thread_store_input (int32_t n_super_cycle, int32_t myID, int32_t nProcs);
		int32_t thread_store_output (int32_t n_super_cycle, int32_t myID, int32_t nProcs);

		void InitHostStore(int32_t n);
		void FreeHostStore(int32_t n);

        void CudaDummy();

	private:
		int64_t * * pointer_list;
		int32_t n_pointer;
		int32_t store_size;
		int32_t size_device;
		int32_t size_debug;
		int32_t size_develop_a;
		int32_t n_store_in;
		int32_t n_store_worker;
		int32_t n_store_host;
		size_t size_temp;

		int32_t n_store_out;

		volatile mem_info * stat_mem_in;
		volatile mem_info * stat_mem_out;
		volatile mem_info * stat_mem_worker;
		volatile mem_info * stat_mem_store;

		CUdeviceptr * d_store;
		CUdeviceptr * d_in;
		CUdeviceptr * d_out;
		CUdeviceptr * d_worker;

		double * h_store [n_h_store_max];

		// algorithm specific data
		CUdeviceptr d_debug;
		CUdeviceptr d_sum;
		CUdeviceptr d_temp;
		CUdeviceptr d_prefix_sum;
		CUdeviceptr d_mol_list;

		CUdeviceptr d_sum_u;
		CUdeviceptr d_sum_v;
		CUdeviceptr d_sum_vel;

		CUdeviceptr d_part_u;
		CUdeviceptr d_part_v;
		CUdeviceptr d_part_vel;
		CUdeviceptr d_part_r_u;
		CUdeviceptr d_part_r_v;

		CUdeviceptr d_result;
		CUdeviceptr d_sum_N;
		CUdeviceptr d_tmp_f;

		double * h_sum_u;
		double * h_sum_v;
		double * h_sum_vel;
		int64_t * h_sum_N;
		int32_t n_cycle_sample;


		CUdeviceptr d_develop_a;

		// end of algorithm specific data

		cudaStream_t stream_in;
		cudaStream_t stream_out;
		cudaStream_t stream_worker;

		cudaEvent_t worker_event [n_worker_event];

		int32_t worker_threads_per_block;
		int32_t worker_n_block;

		int32_t i_gpu;
		int32_t n_worker;
		int32_t n_worker_total;
		int32_t n_part;
		int32_t order_in;
		int32_t order_out;

		MPI_Request request_results[2];
		bool outstanding_results;
		int32_t outstanding_results_n_cycle_sample;
		double vdata_send_a[3*1024];
		double vdata_rec_a[3*1024];
		int64_t vdata_send_b[1024];
		int64_t vdata_rec_b[1024];
		int32_t my_id;
		int32_t n_procs;
		int32_t n_rails;

		char * MyNewCatchChar(const char * file, int32_t line, int64_t size);

		// control logic
		int32_t part_in_present_wait (int32_t i_part, int32_t i_cycle);
		void part_out_ready_wait (int32_t i_part, int32_t i_center);
		void update_mem_info(int32_t *islot, volatile mem_info *mem, int32_t * i_event);
		
		// block management
		int64_t get_block_size_host(double * p_data);
		int64_t get_block_size_device(double * p_data);
        int32_t block_check(double *p_data,int32_t pos);
		int32_t block_check_device (double* data, int32_t pos);
        int32_t block_get_nm(double *p_data);
		int64_t get_block_nm_device(double * p_data);
		int32_t block_mark(double * p_data);
		int32_t block_markcheck(double * p_data);

        // I/O routines
		int32_t FileFieldsToMem(int64_t *dat, char *FileName, int32_t iMax, int64_t *size);
		int64_t MemToFile(int64_t * dat, int64_t n, char * FileName, int32_t newfile);

		// cuda routines
		void cudaCheckError(int32_t line, const char *file);
        int32_t InitCuda();
        void FreeCuda();
		void caller_worker(double **in, double **out, int32_t i_part, int32_t i_super_cycle, int32_t order_in, int32_t order_out, int32_t iworker, int32_t nworker, cudaStream_t *stream, int32_t gridSize, int32_t blockSize, int32_t myID);
};
