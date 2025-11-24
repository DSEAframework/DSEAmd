// Data Streaming for Explicit Algorithms - DSEA

#include <dsea.h>

using namespace std;

// this routine waits for output slots to become free and claims the slots for a block
// todo: move into class
void DS::part_out_ready_wait (int32_t i_part, int32_t i_center) {
	for (int32_t j=-order_out;j<order_out+1;j++) {
		int32_t k=i_part+j;
		if ((k>=0) && (k<=n_part-1)) {
			int32_t islot=i_center+j;	// memory slot
			if (islot<0) islot+=n_store_out;
			if (islot>=n_store_out) islot-=n_store_out;

			if (stat_mem_out[islot].state==mem_state_free) {
				// no waiting required
				stat_mem_out[islot].i_part=k;
			}
			else if (stat_mem_out[islot].state==mem_state_bussy) {
				if (i_part==0) {
				// at the start of the super cycle, the first part needs all memory slots in state free
					// int nwait=0;
					while (stat_mem_out[islot].state!=mem_state_free) {
						}	// wait till slot is free

					stat_mem_out[islot].i_part=k;
				}
				if (j==order_out) {//&&(i_part<n_part-order_out+1)) {
					int nwait=0;
					while (stat_mem_out[islot].state!=mem_state_free) {
						}	// wait till slot is free

					stat_mem_out[islot].i_part=k;
				}
			}
		}
	}
}


// this routine waits till defined part and required neighbours are present in the input buffers
// returns location of part
// todo: move into class
int32_t DS::part_in_present_wait (int32_t i_part, int32_t i_cycle) {
	//bool data_present=true;
	int32_t ifirst=-1;
	int32_t last_test=-1;
	for (int32_t j=-order_in;j<order_in+1;j++) {
		int32_t k=i_part+j;
//		cout << "part_present " << k << endl;
		bool k_found=false;
		while (k_found==false) {
			if ((k>=0) && (k<=n_part-1)) {
				for (int32_t i=0;i<n_store_in;i++) {
					last_test=i;
					if ((stat_mem_in[i].i_part==k)&&(stat_mem_in[i].i_cycle==i_cycle)&&((stat_mem_in[i].state==mem_state_ready)||(stat_mem_in[i].state==mem_state_ready_b)||(stat_mem_in[i].state==mem_state_bussy))) {
						k_found=true;
						break;
					}
				}
			}
			else {
				k_found=true;
			}
		}
		if ((k>=0) && (k<=n_part-1) && (k_found==true) && (ifirst==-1) && (k==i_part)) {
			ifirst=last_test;
		}
	}
	// cout << endl;
	return ifirst;
}


void DS::update_mem_info (int32_t * islot, volatile mem_info * mem, int32_t * i_event) {
	if (mem[0].type==mem_type_input) {
		for (int i=0;i<16;i++) {
			int32_t j=islot[i];
			if (j!=-1) {
				mem[j].state=mem_state_bussy; // todo: not sure if this is necessary!
				mem[j].n_use++;
				// check if this is the last usage of this block, if yes add event to wait for
				int32_t max_use=1+2*order_in;

				// parts at the start and end are used less often
				if (mem[j].i_part<order_in) {
					// first parts
					max_use-=order_in-mem[j].i_part;
					// cout << "ia_" << i << " " << j << " " << mem[j].i_part << " " << order_in << endl;
				}
				else if (n_part-mem[j].i_part-1<order_in) {
					// last parts
					max_use-=n_part-mem[j].i_part;
					// cout << "ib_" << i << " " << j << " " << mem[j].i_part << " " << order_in << endl;
				}

				if (mem[j].n_use==max_use) {
					// cout << "input_last_use_" << j << "_" << max_use << "_" << *i_event << endl;
					cudaError_t cer=cudaEventRecord(worker_event[*i_event],stream_worker);			cudaCheckError(__LINE__,__FILE__);
					if (cer==cudaSuccess) {
						mem[j].i_event=*i_event;
						mem[j].i_cycle=-1;	// this might make the check for i_cycle in part_in_present_wait obsolete
						mem[j].i_part=-1;
						(*i_event)++;
					}
				}
			}
		}
	}
	else if (mem[0].type==mem_type_output) {
		for (int i=0;i<16;i++) {
			int32_t j=islot[i];
			if (j!=-1) {
				mem[j].state=mem_state_bussy; // todo: not sure if this is necessary!
				mem[j].n_use++;
				// check if this is the last usage of this block, if yes add event to wait for
				int32_t max_use=1+2*order_out;

				// parts at the start and end are used less often
				if (mem[j].i_part<order_out) {
					// first parts
					max_use-=order_out-mem[j].i_part;
					// cout << "oa_" << i << " " << j << " " << mem[j].i_part << " " << order_out << endl;
				}
				else if (n_part-mem[j].i_part-1<order_out) {
					// last parts
					max_use-=n_part-mem[j].i_part;
					// cout << "ob_" << i << " " << j << " " << mem[j].i_part << " " << order_out << endl;
				}

				if (mem[j].n_use==max_use) {
					// cout << "output_last_use_" << j << "_" << max_use << "_" << *i_event << endl;
					cudaError_t cer=cudaEventRecord(worker_event[*i_event],stream_worker);			cudaCheckError(__LINE__,__FILE__);
					if (cer==cudaSuccess) {
						mem[j].i_event=*i_event;
						mem[j].i_cycle=-1;	// this might make the check for i_cycle in part_in_present_wait obsolete
						mem[j].i_part=-1;
						(*i_event)++;
					}
				}
			}
		}
	}
	else {
		cout << "not_impl_update_mem_info" << endl;
	}
}

int32_t DS::thread_main (
				int32_t n_super_cycle,
				int32_t order_in,
				int32_t order_out,
				int32_t myID) {

	// current part of worker
	int32_t * w_part = new int32_t [n_worker];
	int32_t * w_i_in = new int32_t [n_worker];
	int32_t * w_i_out = new int32_t [n_worker];

	int32_t last_worker_mem_out_center=0;

	int32_t i_worker_event=0;

	for (int32_t i_super_cycle=0;i_super_cycle<n_super_cycle;i_super_cycle++) {
		if (myID==0) cout << "super_cycle: " << i_super_cycle << endl;

		for (int32_t i=0;i<n_worker;i++) {
				w_part[i]=0;
				w_i_in[i]=0;
				w_i_out[i]=0;
		}

		int32_t n_stage=n_part;//*n_worker;
		bool any_workers_active=true;

		int32_t i_stage=0;
		while (any_workers_active==true) {
			any_workers_active=false;
			// cout << " i_stage " << i_stage <<  endl;
			i_stage++;

			for (int32_t i_worker=0;i_worker<n_worker;i_worker++) {
				int32_t part_to_process=w_part[i_worker];
				// current worker wants to process part part_tp_process

				bool worker_first = (i_worker==0);
				bool worker_last = (i_worker==(n_worker-1));

				// determine input arrays

				double * p_in [16];
				double * p_out [16];
				int32_t islot_in [16];
				int32_t islot_out [16];

				for (int32_t i=0;i<16;i++) {
					p_in[i]={(double*)-1};
					p_out[i]={(double*)-1};
					islot_in[i]=-1;
					islot_out[i]=-1;
				}

				bool worker_active=true;

				if (w_part[i_worker]==my_n_part) {
					worker_active=false;	// worker done
				}
				else {
					if (i_worker>0) {
						if ((w_part[i_worker-1]==my_n_part)||(w_part[i_worker-1]-part_to_process>order_out+1)) {
							// worker active
						}
						else {
							worker_active=false;
						}
					}
				}

				// cout << "i_worker: " << i_worker << " part: " << part_to_process << " active: " << worker_active << endl;

				if (worker_active==true) {
					any_workers_active=true;
					// determine input arrays
					if (worker_first) {
						// first worker uses input from d_in
						int res=part_in_present_wait(part_to_process,i_super_cycle);

						int32_t ibuf=0;
						int32_t slot=res;
						// cout << "slot_to_in_buf:"<< slot << "_" << ibuf << "_" << k <<endl;
						islot_in[ibuf]=slot;
						p_in[ibuf]=(double*)d_in[slot];
						ibuf++;

						for (int32_t j=1;j<order_in+1;j++) {
							int32_t k=part_to_process-j;	// index of input block
							if ((k>=0)&&(k<=n_part-1)) {
								slot=res-j;
								if (slot<0) slot+=n_store_in;
								// cout << "slot_to_in_buf:"<< slot << "_" << ibuf << endl;
								islot_in[ibuf]=slot;
								p_in[ibuf]=(double*)d_in[slot];
							}
							ibuf++;
							k=part_to_process+j;	// index of input block
							if ((k>=0)&&(k<=n_part-1)) {
								slot=res+j;
								if (slot>=n_store_in) slot-=n_store_in;
								// cout << "slot_to_in_buf:"<< slot << "_" << ibuf << "_" << k <<endl;
								islot_in[ibuf]=slot;
								p_in[ibuf]=(double*)d_in[slot];
							}
							ibuf++;
						}
					}
					else {
						int32_t ibuf=0;
						int32_t worker_mem_in_center=w_i_in[i_worker];

						// center part
						int32_t slot=worker_mem_in_center;
						islot_in[ibuf]=slot;
						p_in[ibuf]=(double*)d_worker[(i_worker-1)*n_store_worker+slot];
						ibuf++;

						for (int32_t j=1;j<order_in+1;j++) {
							int32_t k=part_to_process-j;
							if ((k>=0)&&(k<=n_part-1)) {
								slot=worker_mem_in_center-j;
								if (slot<0) slot+=n_store_worker;
								// cout << "slot_to_out_buf:"<< slot << "_" << ibuf << endl;
								islot_in[ibuf]=slot;
								p_in[ibuf]=(double*)d_worker[(i_worker-1)*n_store_worker+slot];
							}
							ibuf++;
							k=part_to_process+j;
							if ((k>=0)&&(k<=n_part-1)) {
								slot=worker_mem_in_center+j;
								// if (slot<0) slot+=n_store_out;
								if (slot>=n_store_worker) slot-=n_store_worker;
								// cout << "slot_to_out_buf:"<< slot << "_" << ibuf << endl;
								islot_in[ibuf]=slot;
								p_in[ibuf]=(double*)d_worker[(i_worker-1)*n_store_worker+slot];
							}
							ibuf++;
						}
						worker_mem_in_center++;
						if (worker_mem_in_center==n_store_worker) worker_mem_in_center=0;
						w_i_in[i_worker]=worker_mem_in_center;
					}

					// determine output arrays
					if (worker_last) {
						part_out_ready_wait(part_to_process,last_worker_mem_out_center);
						int32_t ibuf=0;
						int32_t slot=last_worker_mem_out_center;
						// cout << "slot_to_out_buf:"<< slot << "_" << ibuf << "_" << k <<endl;
						islot_out[ibuf]=slot;
						p_out[ibuf]=(double*)d_out[slot];
						ibuf++;

						for (int32_t j=1;j<order_out+1;j++) {
							int32_t k=part_to_process-j;
							if ((k>=0)&&(k<=n_part-1)) {
								slot=last_worker_mem_out_center-j;
								if (slot<0) slot+=n_store_out;
								// cout << "slot_to_out_buf:"<< slot << "_" << ibuf << endl;
								islot_out[ibuf]=slot;
								p_out[ibuf]=(double*)d_out[slot];
							}
							ibuf++;
							k=part_to_process+j;
							if ((k>=0)&&(k<=n_part-1)) {
								slot=last_worker_mem_out_center+j;
								if (slot>=n_store_out) slot-=n_store_out;
								// cout << "slot_to_out_buf:"<< slot << "_" << ibuf << endl;
								islot_out[ibuf]=slot;
								p_out[ibuf]=(double*)d_out[slot];
							}
							ibuf++;
						}
						last_worker_mem_out_center++;
						if (last_worker_mem_out_center==n_store_out) last_worker_mem_out_center=0;
					}
					else {
						int32_t ibuf=0;
						int32_t worker_mem_out_center=w_i_out[i_worker];

						// center part
						int32_t slot=worker_mem_out_center;
						islot_out[ibuf]=slot;
						p_out[ibuf]=(double*)d_worker[i_worker*n_store_worker+slot];
						ibuf++;

						for (int32_t j=1;j<order_out+1;j++) {
							int32_t k=part_to_process-j;
							if ((k>=0)&&(k<=n_part-1)) {
								slot=worker_mem_out_center-j;
								if (slot<0) slot+=n_store_worker;
								// cout << "slot_to_out_buf:"<< slot << "_" << ibuf << endl;
								islot_out[ibuf]=slot;
								p_out[ibuf]=(double*)d_worker[i_worker*n_store_worker+slot];
							}
							ibuf++;
							k=part_to_process+j;
							if ((k>=0)&&(k<=n_part-1)) {
								slot=worker_mem_out_center+j;
								// if (slot<0) slot+=n_store_out;
								if (slot>=n_store_worker) slot-=n_store_worker;
								// cout << "slot_to_out_buf:"<< slot << "_" << ibuf << endl;
								islot_out[ibuf]=slot;
								p_out[ibuf]=(double*)d_worker[i_worker*n_store_worker+slot];
							}
							ibuf++;
						}
						worker_mem_out_center++;
						if (worker_mem_out_center==n_store_worker) worker_mem_out_center=0;
						w_i_out[i_worker]=worker_mem_out_center;
					}

					// enqueue kernel
					// for (int32_t i=0;i<8;i++) {
						// cout << p_in[i] << " " << p_out[i] << endl;
						// cout << i << "_" << islot_in[i] << "_" << islot_out[i] << endl;
					// }

					// reference version - no optimizations
					ref_caller_worker (p_in,p_out,part_to_process,i_super_cycle,order_in,order_out,i_worker,n_worker,&stream_worker,worker_threads_per_block,worker_n_block,myID);

					// optimized version - DSEA paper
					// b_caller_worker (p_in,p_out,part_to_process,i_super_cycle,order_in,order_out,i_worker,n_worker,&stream_worker,worker_threads_per_block,worker_n_block,myID);

					// visualization
					if (false)
					if (my_id==0) {
					if (worker_first) {
					// if (i_super_cycle > 0)
					if (i_super_cycle % 10 == 0) {
						// cout << "output_pre_" << part_to_process << endl;
						// caller_output_vtk(p_in[0],(double*)d_visual,&stream_worker,worker_threads_per_block,worker_n_block,myID,i_super_cycle);
						caller_output_vtk_rectilinear(p_in[0],(double*)d_visual,&stream_worker,worker_threads_per_block,worker_n_block,myID,i_super_cycle,part_to_process);
						// cout << "output_post" << endl;
					}
					}
					}

					// record event
					if (worker_first) {
						update_mem_info(islot_in,stat_mem_in,&i_worker_event);
						if (i_worker_event==n_worker_event) i_worker_event=0;
					}
					if (worker_last) {
						update_mem_info(islot_out,stat_mem_out,&i_worker_event);
						if (i_worker_event==n_worker_event) i_worker_event=0;
					}

		
					for (int i=0;i<n_store_in;i++) {
						// cout << i << "_use_in_out_"<< "_" << stat_mem_in[i].n_use << "_" << stat_mem_out[i].n_use << "  part_in_out_" << stat_mem_in[i].i_part << "_" << stat_mem_out[i].i_part <<endl;
					}

					w_part[i_worker]++;
				}
			}
		}
	}

	return 0;
}