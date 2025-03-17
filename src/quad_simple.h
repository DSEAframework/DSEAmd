#ifndef UCX_MR_SIMPLE_H
#define UCX_MR_SIMPLE_H

#include <cuda_runtime.h>
#include <cuda.h>
#include <driver_types.h>

#include <math.h>

#include <ucp/api/ucp.h>

#include "ucx_mr_types.h"
#include "ucx_mr_wait.h"

#include <iostream>
using namespace std;


void calculate_msg_sizes_quad(size_t msg_size, float split_ratio, unsigned int element_size, size_t *msg_sizes)
{
  size_t nof_e, nof_e_stay, nof_e_split, nof_e_rest;

  nof_e = (size_t)roundf((msg_size / element_size));

  DEBUG_PRINT("Full message size: %ld Bytes, nof_e: %ld\n", msg_size, nof_e);

  nof_e_stay = (size_t)roundf((msg_size / element_size) * (1 - split_ratio));
  nof_e_split = (nof_e - nof_e_stay) / 3;
  nof_e_rest = (nof_e - nof_e_stay) % 3;

  msg_sizes[0] = nof_e_stay * element_size;
  msg_sizes[1] = nof_e_split * element_size;
  msg_sizes[2] = nof_e_split * element_size;
  msg_sizes[3] = (nof_e_split + nof_e_rest) * element_size;

  DEBUG_PRINT(
      "Calculated: \n msg_sizes0: %ld Bytes, \n msg_sizes1: %ld Bytes, \n msg_sizes2: %ld Bytes, \n msg_sizes3: %ld Bytes, \n Sum: %ld Bytes\n",
      msg_sizes[0], msg_sizes[1], msg_sizes[2], msg_sizes[3], msg_sizes[0] + msg_sizes[1] + msg_sizes[2] + msg_sizes[3]);
}

void calculate_msg_sizes_tripple(size_t msg_size, float split_ratio, unsigned int element_size, size_t *msg_sizes)
{
  size_t nof_e, nof_e_stay, nof_e_split, nof_e_rest;

  nof_e = (size_t)roundf((msg_size / element_size));

  DEBUG_PRINT("Full message size: %ld Bytes, nof_e: %ld\n", msg_size, nof_e);

  nof_e_stay = (size_t)roundf((msg_size / element_size) * (1 - split_ratio));
  nof_e_split = (nof_e - nof_e_stay) / 2;
  nof_e_rest = (nof_e - nof_e_stay) % 2;

  msg_sizes[0] = nof_e_stay * element_size;
  msg_sizes[1] = nof_e_split * element_size;
  msg_sizes[2] = (nof_e_split + nof_e_rest) * element_size;

  DEBUG_PRINT(
      "Calculated: \n msg_sizes0: %ld Bytes, \n msg_sizes1: %ld Bytes, \n msg_sizes2: %ld Bytes, \n msg_sizes3: %ld Bytes, \n Sum: %ld Bytes\n",
      msg_sizes[0], msg_sizes[1], msg_sizes[2], msg_sizes[0] + msg_sizes[1] + msg_sizes[2]);
}

const uint32_t crc32_tab[] = {
	0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419, 0x706af48f,
	0xe963a535, 0x9e6495a3,	0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988,
	0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91, 0x1db71064, 0x6ab020f2,
	0xf3b97148, 0x84be41de,	0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
	0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec,	0x14015c4f, 0x63066cd9,
	0xfa0f3d63, 0x8d080df5,	0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
	0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b,	0x35b5a8fa, 0x42b2986c,
	0xdbbbc9d6, 0xacbcf940,	0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
	0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423,
	0xcfba9599, 0xb8bda50f, 0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
	0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d,	0x76dc4190, 0x01db7106,
	0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
	0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818, 0x7f6a0dbb, 0x086d3d2d,
	0x91646c97, 0xe6635c01, 0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e,
	0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457, 0x65b0d9c6, 0x12b7e950,
	0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
	0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2, 0x4adfa541, 0x3dd895d7,
	0xa4d1c46d, 0xd3d6f4fb, 0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
	0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9, 0x5005713c, 0x270241aa,
	0xbe0b1010, 0xc90c2086, 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
	0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17, 0x2eb40d81,
	0xb7bd5c3b, 0xc0ba6cad, 0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a,
	0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683, 0xe3630b12, 0x94643b84,
	0x0d6d6a3e, 0x7a6a5aa8, 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
	0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb,
	0x196c3671, 0x6e6b06e7, 0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc,
	0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5, 0xd6d6a3e8, 0xa1d1937e,
	0x38d8c2c4, 0x4fdff252, 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
	0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60, 0xdf60efc3, 0xa867df55,
	0x316e8eef, 0x4669be79, 0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
	0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f, 0xc5ba3bbe, 0xb2bd0b28,
	0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
	0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a, 0x9c0906a9, 0xeb0e363f,
	0x72076785, 0x05005713, 0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38,
	0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21, 0x86d3d2d4, 0xf1d4e242,
	0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
	0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c, 0x8f659eff, 0xf862ae69,
	0x616bffd3, 0x166ccf45, 0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2,
	0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db, 0xaed16a4a, 0xd9d65adc,
	0x40df0b66, 0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
	0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6, 0xbad03605, 0xcdd70693,
	0x54de5729, 0x23d967bf, 0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
	0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d
};

uint32_t crc32(const void *buf, size_t size)
{
    const uint8_t *p = (const uint8_t * )buf;
    uint32_t crc;

    crc = ~0U;
    while (size--)
    crc = crc32_tab[(crc ^ *p++) & 0xFF] ^ (crc >> 8);
    return crc ^ ~0U;
}

ucs_status_t
ucx_mr_quad_split_send_simple(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
                       void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                       void *buffer1, ucs_memory_type_t mem_type1, int device1,
                       void *buffer2, ucs_memory_type_t mem_type2, int device2,
                       void *buffer3, ucs_memory_type_t mem_type3, int device3)
{

// // determine crc of data
// char * p_debug=new char [msg_size0];
// cudaMemcpy((void*)p_debug,(const void*)buffer0,msg_size0,cudaMemcpyDeviceToHost);//		cudaCheckError(__LINE__,__FILE__);
// cout << "crc_send:_" << msg_size0 << "_" << crc32(p_debug,msg_size0) << endl << flush;
// delete [] p_debug;

  // split_ratio : factor how much is sent to peer (0.39 -> 39% is sent to peer devices)
  // element_size: nof bytes of splitted memory needs to be multiple of element_size
  // buf_size1   : maximum possible splittable memory
  cudaError_t c_status;
  ucs_status_t status[4] = {UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE};

  struct ucx_context *request[4] = {NULL,NULL,NULL,NULL};
  ucp_request_param_t send_param[4];

  size_t msg_sizes[4];

  int devices[4]                     = {device0, device1, device2, device3};
  ucs_memory_type_t mem_types[4]     = {mem_type0, mem_type1, mem_type2, mem_type3};
  void *buffer[4]                    = {buffer0, buffer1, buffer2, buffer3};

  ucp_ep_h *ep              = mr_ctx->ep;
  ucp_worker_h *worker      = mr_ctx->worker;

//   cudaStream_t *cuda_streams  = mr_ctx->cuda_streams;
//   cudaEvent_t *cuda_events    = mr_ctx->cuda_events;

  for (int i = 0; i < 4; ++i)
  {
    ucx_mr_write_send_param(&send_param[i], mem_types[i]);
  }

  calculate_msg_sizes_quad(msg_size0, split_ratio, element_size, msg_sizes);
  //calculate_msg_sizes_quad_alt(msg_size0, split_ratio, element_size, msg_sizes);

  cudaSetDevice(devices[0]);
  request[0] = (ucx_context*)ucp_tag_send_nbx(ep[0], buffer[0], msg_sizes[0], tag, &send_param[0]);

  
  size_t offset = 0;
  for (int i = 0; i < 3; ++i)
  {

    offset += msg_sizes[i];

    DEBUG_PRINT("Offset: %ld\n", offset);
    cudaSetDevice(devices[0]);
    c_status = cudaMemcpyPeer(buffer[i + 1], devices[i + 1], buffer0 + offset, device0, msg_sizes[i + 1]);
    if (c_status != cudaSuccess)
    {
      DEBUG_PRINT("CudaMemcpy Peer Error!\n");
      return UCS_ERR_REJECTED;
    }
    cudaSetDevice(devices[i+1]);
    // cout << "send_state_k_" << i << endl;
    request[i+1] = (ucx_context*)ucp_tag_send_nbx(ep[i+1], buffer[i+1], msg_sizes[i+1], tag, &send_param[i+1]);
  }


//   for (int i = 0; i < 4; ++i)
//   {
// // cout << "send_state_a_" << i << endl;

//     cudaSetDevice(devices[i]);
//     status[i] = ucx_wait(worker[i], request[i], "receive", data_msg_str);


// // cout << "send_state_b_" << i << "_" << status[i] << endl;

//     if (status[i] != UCS_OK) {
//     //   ERROR_PRINT("ucx send error!\n");
//       return UCS_ERR_REJECTED;
//     }
//   }


	int n_ok=0;
	while (n_ok!=4) {
		for (int i = 0; i < 4; ++i) {
			if (status[i]!=UCS_OK) {
				cudaSetDevice(devices[i]);
				status[i] = ucx_wait_max(worker[i], (ucx_context*)request[i], "receive", data_msg_str,10);
				if (status[i] != UCS_OK) {
					//   ERROR_PRINT("ucx send error!\n");
					// return UCS_ERR_REJECTED;
					}
					else {
						n_ok++;
				}
		  }
		}
	}



    cudaSetDevice(devices[0]);
  
  return UCS_OK;
}


ucs_status_t
ucx_mr_tripple_split_send_simple(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
                       void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                       void *buffer1, ucs_memory_type_t mem_type1, int device1,
                       void *buffer2, ucs_memory_type_t mem_type2, int device2)
{

// // determine crc of data
// char * p_debug=new char [msg_size0];
// cudaMemcpy((void*)p_debug,(const void*)buffer0,msg_size0,cudaMemcpyDeviceToHost);//		cudaCheckError(__LINE__,__FILE__);
// cout << "crc_send:_" << msg_size0 << "_" << crc32(p_debug,msg_size0) << endl << flush;
// delete [] p_debug;


  // split_ratio : factor how much is sent to peer (0.39 -> 39% is sent to peer devices)
  // element_size: nof bytes of splitted memory needs to be multiple of element_size
  // buf_size1   : maximum possible splittable memory
  cudaError_t c_status;
  ucs_status_t status[3] = {UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE};

  struct ucx_context *request[3] = {NULL};
  ucp_request_param_t send_param[3];

  size_t msg_sizes[3];

  int devices[3]                     = {device0, device1, device2};
  ucs_memory_type_t mem_types[3]     = {mem_type0, mem_type1, mem_type2};
  void *buffer[3]                    = {buffer0, buffer1, buffer2};

  ucp_ep_h *ep              = mr_ctx->ep;
  ucp_worker_h *worker      = mr_ctx->worker;

  for (int i = 0; i < 3; ++i)
  {
    ucx_mr_write_send_param(&send_param[i], mem_types[i]);
  }

  calculate_msg_sizes_tripple(msg_size0, split_ratio, element_size, msg_sizes);
  
  cudaSetDevice(devices[0]);
  request[0] = (ucx_context*)ucp_tag_send_nbx(ep[0], buffer[0], msg_sizes[0], tag, &send_param[0]);

  
  size_t offset = 0;
  for (int i = 0; i < 2; ++i)
  {

    offset += msg_sizes[i];

    DEBUG_PRINT("Offset: %ld\n", offset);
    cudaSetDevice(devices[0]);
    c_status = cudaMemcpyPeer(buffer[i + 1], devices[i + 1], buffer0 + offset, device0, msg_sizes[i + 1]);
    if (c_status != cudaSuccess)
    {
      cout << "prioblem_cudaMemcpyPeer" << endl;


	// cudaDeviceSynchronize();
	// cudaError_t e = cudaGetLastError();
	// if (e != cudaSuccess) {
	// 		printf("Cuda failure %s:%d: '%s'\n", __FILE__,__LINE__,cudaGetErrorString(e));
	// 		exit(EXIT_FAILURE);
	// }


      DEBUG_PRINT("CudaMemcpy Peer Error!\n");
      return UCS_ERR_REJECTED;
    }
    cudaSetDevice(devices[i+1]);
    // cout << "send_state_k_" << i << endl;
    request[i+1] = (ucx_context*)ucp_tag_send_nbx(ep[i+1], buffer[i+1], msg_sizes[i+1], tag+i+1, &send_param[i+1]);
  }


//   for (int i = 0; i < 3; ++i)
//   {
// // cout << "send_state_a_" << i << endl;

//     cudaSetDevice(devices[i]);
//     status[i] = ucx_wait(worker[i], request[i], "receive", data_msg_str);
//     cout << status[i] << " " << request[i] << endl;
// // cout << "send_state_b_" << i << "_" << status[i] << endl;

//     if (status[i] != UCS_OK) {
//     //   ERROR_PRINT("ucx send error!\n");
//       return UCS_ERR_REJECTED;
//     }
//   }




	int n_ok=0;
	while (n_ok!=3) {
		for (int i = 0; i < 3; ++i) {
			if (status[i]!=UCS_OK) {
				cudaSetDevice(devices[i]);
				status[i] = ucx_wait_max(worker[i], (ucx_context*)request[i], "receive", data_msg_str,10);
				if (status[i] != UCS_OK) {
					//   ERROR_PRINT("ucx send error!\n");
					// return UCS_ERR_REJECTED;
					}
					else {
						n_ok++;
				}
		  }
		}
	}

    cudaSetDevice(devices[0]);
  
  return UCS_OK;
}


ucs_status_t
ucx_mr_quad_split_recv_simple(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
                       void *buffer0, size_t buffer_size, ucs_memory_type_t mem_type0, int device0,
                       void *buffer1, ucs_memory_type_t mem_type1, int device1,
                       void *buffer2, ucs_memory_type_t mem_type2, int device2,
                       void *buffer3, ucs_memory_type_t mem_type3, int device3)
{
    // cout << "enter_rec" << endl;
  // 0 is the collecting device
  cudaError_t c_status;
  ucs_status_t status[4] = {UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE};

  struct ucx_context *request[4] = {NULL, NULL, NULL, NULL};
  ucp_request_param_t recv_param[4];

  ucp_tag_recv_info_t info_tag[4];
  ucp_tag_message_h msg_tag[4] = {NULL, NULL, NULL, NULL};
  size_t msg_sizes[4] = {0};


  int devices[4]                  = {device0, device1, device2, device3};
  ucs_memory_type_t mem_types[4]  = {mem_type0, mem_type1, mem_type2, mem_type3};
  void *buffer[4]                 = {buffer0, buffer1, buffer2, buffer3};
  size_t buffer_sizes[4]          = {buffer_size, buffer_size, buffer_size, buffer_size};

  ucp_worker_h *worker = mr_ctx->worker;

//   cudaStream_t *cuda_streams = mr_ctx->cuda_streams;
  //cudaEvent_t *cuda_events = mr_ctx->cuda_events;

  for (int i = 0; i < 4; ++i)
  {
    ucx_mr_write_recv_param(&recv_param[i], mem_types[i]);
  }
// cout << "rec_state_a" << endl;

  for (;;)
  {
    for (int i = 0; i < 4; ++i)
    {
      if (msg_tag[i] == NULL)
      {
        cudaSetDevice(devices[i]);
        msg_tag[i] = ucp_tag_probe_nb(worker[i], tag, tag_mask, 1, &info_tag[i]);

        if (msg_tag[i] != NULL)
        {
          request[i] = (ucx_context*)ucp_tag_msg_recv_nbx(worker[i], buffer[i], buffer_sizes[i], msg_tag[i], &recv_param[i]);
          msg_sizes[i] = info_tag[i].length;
        //   cout << "rec_state_c" << "_" << i << endl;
        }
        ucp_worker_progress(worker[i]);
      }
    }

    if (msg_tag[0] != NULL && msg_tag[1] != NULL && msg_tag[2] != NULL && msg_tag[3] != NULL)
      break;
  }
//   cout << "rec_state_b" << endl;

  size_t offset = 0;
  for (int i = 1; i < 4; ++i)
  {
    cudaSetDevice(devices[i]);
    status[i] = ucx_wait(worker[i], request[i], "receive", data_msg_str);
    if (status[i] != UCS_OK) {
    //   ERROR_PRINT("ucx send error!\n");
      return UCS_ERR_REJECTED;
    }
    
    offset += msg_sizes[i-1];
        
    c_status = cudaMemcpyPeer((buffer0 + offset), devices[0], buffer[i], devices[i], msg_sizes[i]);
    if (c_status != cudaSuccess)
    {
    //   ERROR_PRINT("cudaMemcpyPeer Failed!\n");
      return UCS_ERR_REJECTED;
    }
  }

//   cout << "rec_state_c" << endl;

  cudaSetDevice(devices[0]);
  status[0] = ucx_wait(worker[0], request[0], "receive", data_msg_str);
  if (status[0] != UCS_OK) {
    // ERROR_PRINT("ucx send error!\n");
    return UCS_ERR_REJECTED;
  }

// cout << "rec_state_d" << endl;

    // size_t total_msg_size=msg_sizes[0]+msg_sizes[1]+msg_sizes[2]+msg_sizes[3];

    // char * p_debug=new char [total_msg_size];

    // cudaMemcpy((void*)p_debug,(const void*)buffer0,total_msg_size,cudaMemcpyDeviceToHost);//		cudaCheckError(__LINE__,__FILE__);

    // cout << "crc_rec:_" << total_msg_size << "_" << crc32(p_debug,total_msg_size) << endl << flush ;
    // delete [] p_debug;


  return UCS_OK;
}


ucs_status_t
ucx_mr_tripple_split_recv_simple(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
                       void *buffer0, size_t buffer_size, ucs_memory_type_t mem_type0, int device0,
                       void *buffer1, ucs_memory_type_t mem_type1, int device1,
                       void *buffer2, ucs_memory_type_t mem_type2, int device2)
{
  cudaError_t c_status;
  ucs_status_t status[3] = {UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE, UCS_ERR_NO_MESSAGE};

  struct ucx_context *request[3] = {NULL, NULL, NULL};
  ucp_request_param_t recv_param[3];

  ucp_tag_recv_info_t info_tag[3];
  ucp_tag_message_h msg_tag[3] = { NULL, NULL, NULL};
  size_t msg_sizes[3] = {0};


  int devices[3]                  = {device0, device1, device2};
  ucs_memory_type_t mem_types[3]  = {mem_type0, mem_type1, mem_type2};
  void *buffer[3]                 = {buffer0, buffer1, buffer2};
  size_t buffer_sizes[3]          = {buffer_size, buffer_size, buffer_size};

  ucp_worker_h *worker = mr_ctx->worker;

//   cudaStream_t *cuda_streams = mr_ctx->cuda_streams;
  //cudaEvent_t *cuda_events = mr_ctx->cuda_events;

  for (int i = 0; i < 3; ++i)
  {
    ucx_mr_write_recv_param(&recv_param[i], mem_types[i]);
  }

  for (;;)
  {
    for (int i = 0; i < 3; ++i)
    {
      if (msg_tag[i] == NULL)
      {
        cudaSetDevice(devices[i]);
        msg_tag[i] = ucp_tag_probe_nb(worker[i], tag+i, tag_mask, 1, &info_tag[i]);

        if (msg_tag[i] != NULL)
        {
          request[i] = (ucx_context*)ucp_tag_msg_recv_nbx(worker[i], buffer[i], buffer_sizes[i], msg_tag[i], &recv_param[i]);
          msg_sizes[i] = info_tag[i].length;
        }
        ucp_worker_progress(worker[i]);
      }
    }

    if (msg_tag[0] != NULL && msg_tag[1] != NULL && msg_tag[2] != NULL )
      break;
  }

  size_t offset = 0;
  for (int i = 1; i < 3; ++i)
  {
    cudaSetDevice(devices[i]);
    status[i] = ucx_wait(worker[i], request[i], "receive", data_msg_str);
    if (status[i] != UCS_OK) {
    //   ERROR_PRINT("ucx send error!\n");
      return UCS_ERR_REJECTED;
    }
    
    offset += msg_sizes[i-1];
        
    c_status = cudaMemcpyPeer((buffer0 + offset), devices[0], buffer[i], devices[i], msg_sizes[i]);
    if (c_status != cudaSuccess)
    {
    //   ERROR_PRINT("cudaMemcpyPeer Failed!\n");
      return UCS_ERR_REJECTED;
    }
  }

  cudaSetDevice(devices[0]);
  status[0] = ucx_wait(worker[0], request[0], "receive", data_msg_str);
  if (status[0] != UCS_OK) {
    // ERROR_PRINT("ucx send error!\n");
    return UCS_ERR_REJECTED;
  }


    // size_t total_msg_size=msg_sizes[0]+msg_sizes[1]+msg_sizes[2]+msg_sizes[3];

    // char * p_debug=new char [total_msg_size];

    // cudaMemcpy((void*)p_debug,(const void*)buffer0,total_msg_size,cudaMemcpyDeviceToHost);//		cudaCheckError(__LINE__,__FILE__);

    // cout << "crc_rec:_" << total_msg_size << "_" << crc32(p_debug,total_msg_size) << endl << flush ;
    // delete [] p_debug;


  return UCS_OK;
}


#endif // UCX_MR_SIMPLE_H