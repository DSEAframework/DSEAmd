// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>
#include <cuda.h>
#include <sys/time.h>

using namespace std;

// constructor
DS::DS (int32_t igpu, int32_t nworker, int32_t npart, int32_t orderin, int32_t orderout, int32_t myid, int32_t nprocs, int32_t nrails) {
	n_worker=nworker;
	i_gpu=igpu;
	n_part=npart;
    order_in=orderin;
    order_out=orderout;
	my_id=myid;
	n_procs=nprocs;
	n_rails=nrails;
	outstanding_results=false;

	n_worker_total=nprocs;

	store_size=-1;
	worker_n_block=-1;
// #if defined case_30
// 	store_size = 8*1024*1024;              // [bytes]
// 	worker_n_block=1024*16;
// #elif defined case_63
// 	store_size = 8*1024*1024;              // [bytes]
// 	worker_n_block=1024*16;
// #elif defined case_136
// 	store_size = 20*1024*1024;              // [bytes]
// 	worker_n_block=1024*16;
// #elif defined case_293
// 	store_size = 80*1024*1024;              // [bytes]
// 	worker_n_block=1024*16;
// #elif defined case_630
// 	store_size = 400*1024*1024;              // [bytes]
// 	worker_n_block=1024*32*4;
// #else
// #endif


#if defined case_30
	store_size = 8*1024*1024;              // [bytes]
	worker_n_block=256;
#elif defined case_63
	store_size = 8*1024*1024;              // [bytes]
	worker_n_block=1024*1;
#elif defined case_136
	store_size = 18*1024*1024;              // [bytes]
	worker_n_block=1024*4;
#elif defined case_293
	store_size = 100*1024*1024;              // [bytes]
	worker_n_block=1024*20;
#elif defined case_630
	store_size = 400*1024*1024;              // [bytes]
	worker_n_block=1024*32*2;
#elif defined case_900
	store_size = 760*1024*1024;              // [bytes]
	worker_n_block=1024*32*5;
#else
#endif


	worker_threads_per_block=32;

	size_device = store_size;

	size_temp=NULL;

	n_pointer=0;
	n_store_in=12;
	n_store_out=12;
	n_store_worker=8;
	n_store_host=0;
	pointer_list = new int64_t * [1024];

	int n_store_sum=(n_store_in+n_store_out+(n_worker-1)*n_store_worker)*nprocs;
	// int n_store_sum=(n_store_in+n_store_out+n_worker*n_store_worker)*nprocs;

	// try to keep the number of storage slots small
	while ((double)n_store_sum<(double)(my_n_part*1.1)) {
		n_store_in++;
		n_store_sum=(n_store_in+n_store_out+(n_worker-1)*n_store_worker)*nprocs;
	}

	if (my_id==0) {
		cout << "n_store_in, n_store_out =" << n_store_in << ", " << n_store_out << " n_store_sum: " << n_store_sum << " my_n_part: " << my_n_part << endl;
	}


	// status of buffers
	// int n_mem_in = ds.n_store_in;//n_store_in;//1+1+2*order_in;
	// int n_mem_out = ds.n_store_out;//1+1+2*order_out;

	stat_mem_in = new mem_info [n_store_in];
	for (int i=0;i<n_store_in;i++) {
		stat_mem_in[i].type=mem_type_input;
		stat_mem_in[i].i_part=-1;
		stat_mem_in[i].state=mem_state_free;
		stat_mem_in[i].n_use=0;
		stat_mem_in[i].i_event=-1;
		stat_mem_in[i].i_cycle=-1;
	}

	stat_mem_out = new mem_info [n_store_out];
	for (int i=0;i<n_store_out;i++) {
		stat_mem_out[i].type=mem_type_output;
		stat_mem_out[i].i_part=-1;
		stat_mem_out[i].state=mem_state_free;
		stat_mem_out[i].n_use=0;
		stat_mem_out[i].i_event=-1;
		stat_mem_out[i].i_cycle=-1;
	}

	if (n_worker>1) {
		stat_mem_worker = new mem_info [n_store_worker*(n_worker-1)];
		for (int i=0;i<n_store_worker;i++) {
			stat_mem_worker[i].type=mem_type_worker;
			stat_mem_worker[i].i_part=-1;
			stat_mem_worker[i].state=mem_state_free;
			stat_mem_worker[i].n_use=0;
			stat_mem_worker[i].i_event=-1;
			stat_mem_worker[i].i_cycle=-1;
		}
	}

	InitCuda();

#if defined case_30
	cout << "case_30" << endl;
#elif defined case_63
	cout << "case_63" << endl;
#elif defined case_136
	cout << "case_136" << endl;
#elif defined case_293
	cout << "case_293" << endl;
#elif defined case_630
	cout << "case_630" << endl;
#elif defined case_900
	cout << "case_900" << endl;
#else
	cout << "case_not_defined" << endl;
#endif

}
 
// destructor
DS::~DS() {
	FreeCuda();

	// free data fields
	for (int32_t i_pointer=0;i_pointer<n_pointer;i_pointer++) {
		if (*pointer_list[i_pointer]!=-1) {
			delete [] (long*)*pointer_list[i_pointer];
		}
	}
	delete [] pointer_list;
}

int64_t DS::MyGetTime() {
	// return time since Epoch in micro seconds
	struct timeval timenow;
	gettimeofday(&timenow, NULL);
	return timenow.tv_sec * 1000000 + timenow.tv_usec;
}
