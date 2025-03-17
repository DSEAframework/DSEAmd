// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include <iostream>
#include <cuda.h>

using namespace std;

void DS::cudaCheckError(int32_t line, const char * file) {
	cudaDeviceSynchronize();
	cudaError_t e = cudaGetLastError();
	if (e != cudaSuccess) {
			printf("Cuda failure %s:%d: '%s'\n", file, line,cudaGetErrorString(e));
			exit(EXIT_FAILURE);
	}
}

//int32_t DS::InitGPUMem () {
//}

void DS::CudaDummy() {
	cudaSetDevice(i_gpu);										cudaCheckError(__LINE__,__FILE__);
}

int32_t DS::InitCuda() {
	int32_t nGPU_available=-1;
	cudaGetDeviceCount(&nGPU_available);						cudaCheckError(__LINE__,__FILE__);
	cout << "nGPU_available: " << nGPU_available << endl;
	cout << "using GPU " << i_gpu << endl;

	cudaSetDevice(i_gpu);										cudaCheckError(__LINE__,__FILE__);


	if (store_size%8!=0) cout << "store_size_must_be_multiple_of_8!" << endl;
	if (size_device%8!=0) cout << "size_device_must_be_multiple_of_8!" << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// host memory - pinned
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//	mem_host_size = 1024*1024*1024;
//	mem_host_size *= 8;

//cudaHostAlloc((void **)&mem_host_i64a,store_size,cudaHostAllocDefault);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// device memory
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int64_t n_alloc=0;

	d_store = new CUdeviceptr [n_worker*3];
	for (int32_t is=0;is<n_worker*3;is++) {
		cudaMalloc((void**)&d_store[is],size_device);					cudaCheckError(__LINE__,__FILE__);
		n_alloc+=size_device;
	}

	d_in = new CUdeviceptr [n_store_in];
	for (int32_t is=0;is<n_store_in;is++) {
		cudaMalloc((void**)&d_in[is],size_device);						cudaCheckError(__LINE__,__FILE__);
		n_alloc+=size_device;
	}

	d_out = new CUdeviceptr [n_store_out];
	for (int32_t is=0;is<n_store_out;is++) {
		cudaMalloc((void**)&d_out[is],size_device);						cudaCheckError(__LINE__,__FILE__);
		n_alloc+=size_device;
	}

	if (n_worker>1) {
		d_worker = new CUdeviceptr [n_store_worker*(n_worker-1)];
		for (int32_t is=0;is<n_store_worker*(n_worker-1);is++) {
			cudaMalloc((void**)&d_worker[is],size_device);					cudaCheckError(__LINE__,__FILE__);
			n_alloc+=size_device;
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// algorithm specific device memory
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	size_debug=1024*1024;
	cudaMalloc((void**)&d_debug,size_debug);									cudaCheckError(__LINE__,__FILE__);
	n_alloc+=size_debug;

	int32_t size_sum=1024*1024*sizeof(double);
	cudaMalloc((void**)&d_sum,size_sum);										cudaCheckError(__LINE__,__FILE__);
	n_alloc+=size_debug;

	cudaMalloc((void**)&d_prefix_sum,block_ncc*sizeof(int32_t));				cudaCheckError(__LINE__,__FILE__);
	n_alloc+=size_sum;

	cudaMalloc((void**)&d_mol_list,block_ncc*c_mol_max*2*sizeof(int32_t));		cudaCheckError(__LINE__,__FILE__);
	n_alloc+=block_ncc*c_mol_max*2*sizeof(int32_t);


	cudaMalloc((void**)&d_sum_u,my_n_part*sizeof(double));					cudaCheckError(__LINE__,__FILE__);
	// n_alloc+=worker_threads_per_block*worker_n_block*sizeof(double);
	cudaMalloc((void**)&d_sum_v,my_n_part*sizeof(double));					cudaCheckError(__LINE__,__FILE__);
	// n_alloc+=worker_threads_per_block*worker_n_block*sizeof(double);
	cudaMalloc((void**)&d_sum_vel,my_n_part*sizeof(double));					cudaCheckError(__LINE__,__FILE__);
	// n_alloc+=worker_threads_per_block*worker_n_block*sizeof(double);

	cudaMalloc((void**)&d_sum_N,my_n_part*sizeof(int64_t));															cudaCheckError(__LINE__,__FILE__);
	n_alloc+=my_n_part*sizeof(int64_t);

	cudaMalloc((void**)&d_result,8*sizeof(double));															cudaCheckError(__LINE__,__FILE__);



	cudaMalloc((void**)&d_part_u,worker_threads_per_block*worker_n_block*sizeof(double));															cudaCheckError(__LINE__,__FILE__);
	n_alloc+=worker_threads_per_block*worker_n_block*sizeof(double);
	cudaMalloc((void**)&d_part_v,worker_threads_per_block*worker_n_block*sizeof(double));															cudaCheckError(__LINE__,__FILE__);
	n_alloc+=worker_threads_per_block*worker_n_block*sizeof(double);
	cudaMalloc((void**)&d_part_vel,worker_threads_per_block*worker_n_block*sizeof(double));															cudaCheckError(__LINE__,__FILE__);
	n_alloc+=worker_threads_per_block*worker_n_block*sizeof(double);

	cudaMalloc((void**)&d_part_r_u,worker_threads_per_block*worker_n_block*sizeof(double));															cudaCheckError(__LINE__,__FILE__);
	n_alloc+=worker_threads_per_block*worker_n_block*sizeof(double);
	cudaMalloc((void**)&d_part_r_v,worker_threads_per_block*worker_n_block*sizeof(double));															cudaCheckError(__LINE__,__FILE__);
	n_alloc+=worker_threads_per_block*worker_n_block*sizeof(double);


	cudaMalloc((void**)&d_tmp_f,3*10000000*sizeof(double));															cudaCheckError(__LINE__,__FILE__);
	n_alloc+=3*10000000*sizeof(double);

	cout << (double)n_alloc/1.0e9 << "_#n_alloc" << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// algorithm specific host memory
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cudaHostAlloc((void**)&h_sum_u,my_n_part*sizeof(double),0);				cudaCheckError(__LINE__,__FILE__);
	cudaHostAlloc((void**)&h_sum_v,my_n_part*sizeof(double),0);				cudaCheckError(__LINE__,__FILE__);
	cudaHostAlloc((void**)&h_sum_vel,my_n_part*sizeof(double),0);				cudaCheckError(__LINE__,__FILE__);
	cudaHostAlloc((void**)&h_sum_N,my_n_part*sizeof(int64_t),0);														cudaCheckError(__LINE__,__FILE__);

	

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// streams
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cudaStreamCreateWithFlags(&stream_worker,cudaStreamNonBlocking);		cudaCheckError(__LINE__,__FILE__);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// events
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int i_event=0;i_event<n_worker_event;i_event++) {
		cudaEventCreate(&worker_event[i_event]);
	}

	return 0;
}

void DS::InitHostStore(int32_t n) {
	n_store_host=n;
	for (int32_t i=0;i<n_h_store_max;i++) {
		h_store[i]=(double*)-1;
	}

	for (int32_t i=0;i<n_store_host;i++) {
		cudaHostAlloc((void **)&h_store[i],store_size,cudaHostAllocDefault);
	}
}

void DS::FreeHostStore(int32_t n) {
	for (int32_t i=0;i<n_store_host;i++) {
		cudaFree((void *)h_store[i]);
	}
}

void DS::FreeCuda() {
	for (int32_t is=0;is<n_worker*3;is++) {
		cudaFree((void*)d_store[is]);				cudaCheckError(__LINE__,__FILE__);
	}
	delete [] d_store;

	for (int32_t is=0;is<n_store_in;is++) {
		cudaFree((void*)d_in[is]);					cudaCheckError(__LINE__,__FILE__);
	}
	delete [] d_in;

	for (int32_t is=0;is<n_store_out;is++) {
		cudaFree((void*)d_out[is]);					cudaCheckError(__LINE__,__FILE__);
	}
	delete [] d_out;

	for (int32_t is=0;is<n_store_worker*(n_worker-1);is++) {
		cudaFree((void*)d_worker[is]);					cudaCheckError(__LINE__,__FILE__);
	}
	delete [] d_worker;

}

char * DS::MyNewCatchChar(const char * file, int32_t line, int64_t size) {
	char * ptmp = new (std::nothrow) char [size];			// allocate new memory
	if ((!ptmp)||(ptmp==NULL)) {
		cout << "allocation of memory failed: " << file << " " << line << " " << size << endl;
		exit(1);
	}
	return ptmp;
}
