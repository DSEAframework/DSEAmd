#ifndef UCX_MR_SEND_H
#define UCX_MR_SEND_H

#include <cuda_runtime.h>
#include <cuda.h>
#include <driver_types.h>

#include <math.h>

#include <ucp/api/ucp.h>

#include "ucx_mr_types.h"
#include "ucx_mr_wait.h"

#include <iostream>
using namespace std;

static void
send_handler(void *request, ucs_status_t status, void *ctx)
{
  struct ucx_context *context = (struct ucx_context *)request;
  const char *str = (const char *)ctx;
  context->completed = 1;
  DEBUG_PRINT("[0x%x] send handler called for \"%s\" with status %d (%s)\n",
              (unsigned int)pthread_self(), str, status,
              ucs_status_string(status));
}

void ucx_mr_write_send_param(ucp_request_param_t *send_param, ucs_memory_type_t mem_type)
{
  send_param->op_attr_mask = UCP_OP_ATTR_FIELD_CALLBACK |
                             UCP_OP_ATTR_FIELD_USER_DATA |
                             UCP_OP_ATTR_FIELD_MEMORY_TYPE;
  send_param->cb.send = send_handler;
  send_param->user_data = (void *)data_msg_str;
  send_param->memory_type = mem_type;
}

ucs_status_t
ucx_mr_parallel_send(ucx_mr_context_t *mr_ctx, ucp_tag_t tag,
                     void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                     void *buffer1, size_t msg_size1, ucs_memory_type_t mem_type1, int device1)
{
  ucs_status_t status0;
  ucs_status_t status1;

  struct ucx_context *request0 = NULL;
  struct ucx_context *request1 = NULL;

  ucp_request_param_t send_param0;
  ucp_request_param_t send_param1;

  ucp_ep_h *ep = mr_ctx->ep;
  ucp_worker_h *worker = mr_ctx->worker;

  ucx_mr_write_send_param(&send_param0, mem_type0);
  ucx_mr_write_send_param(&send_param1, mem_type1);

  cudaSetDevice(device0);
  request0 = (struct ucx_context *)ucp_tag_send_nbx(ep[0], buffer0, msg_size0, tag, &send_param0);
  cudaSetDevice(device1);
  request1 = (struct ucx_context *)ucp_tag_send_nbx(ep[1], buffer1, msg_size1, tag, &send_param1);

  cudaSetDevice(device0);
  status0 = ucx_wait(worker[0], request0, "send", data_msg_str);
  cudaSetDevice(device1);
  status1 = ucx_wait(worker[1], request1, "send", data_msg_str);

  if (status0 != UCS_OK)
    return status0;
  if (status1 != UCS_OK)
    return status1;

  return UCS_OK;
}

ucs_status_t
ucx_mr_split_send(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
                  void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                  void *buffer1, ucs_memory_type_t mem_type1, int device1)
{
  // split_ratio : factor how much is sent to peer (0.39 -> 39% is sent to peer device)
  // element_size: nof bytes of splitted memory needs to be multiple of element_size
  // buf_size1   : maximum possible splittable memory

  cudaError_t c_status;
  ucs_status_t status0;
  ucs_status_t status1;

  struct ucx_context *request0 = NULL;
  struct ucx_context *request1 = NULL;
  ucp_request_param_t send_param0;
  ucp_request_param_t send_param1;

  ucp_ep_h *ep = mr_ctx->ep;
  ucp_worker_h *worker = mr_ctx->worker;

  ucx_mr_write_send_param(&send_param0, mem_type0);
  ucx_mr_write_send_param(&send_param1, mem_type1);

  if (mem_type0 != UCS_MEMORY_TYPE_CUDA && mem_type1 != UCS_MEMORY_TYPE_CUDA)
  {
    printf("Split recv only implemented for CUDA Memory\n");
    return UCS_ERR_NOT_IMPLEMENTED;
  }

  size_t tmp;
  size_t msg_size1;
  size_t nof_elements;

  nof_elements = (size_t)roundf((msg_size0 / element_size) * (1 - split_ratio));

  DEBUG_PRINT("Full message size: %d Bytes, nof_elements: %d\n", msg_size0, nof_elements);

  tmp = nof_elements * element_size;
  msg_size1 = msg_size0 - tmp;
  msg_size0 = tmp;  

  int canAccess = 0;

  DEBUG_PRINT("Calculated msg_size0: %d Bytes, msg_size_1: %d, Sum: %d\n", msg_size0, msg_size1, msg_size0 + msg_size1);

  DEBUG_PRINT("cudaMemcpyPeer: dst: %lld, dstDevice: %d, src: %lld, srcDevice: %d\n", buffer1, device1, buffer0 + msg_size0, device0);


  cudaSetDevice(0);
  request0 = (struct ucx_context *)ucp_tag_send_nbx(ep[0], buffer0, msg_size0, tag, &send_param0);
  
  c_status = cudaMemcpyPeer(buffer1, device1, (buffer0 + msg_size0), device0, msg_size1);
  if (c_status != cudaSuccess)
    return UCS_ERR_REJECTED;

  DEBUG_PRINT("CudaMemcpy Peer finished!\n");

  cudaSetDevice(1);
  request1 = (struct ucx_context *)ucp_tag_send_nbx(ep[1], buffer1, msg_size1, tag, &send_param1);


  status0 = ucx_wait(worker[0], request0, "send", data_msg_str);
  status1 = ucx_wait(worker[1], request1, "send", data_msg_str);

  if (status0 != UCS_OK)
    return status0;
  if (status1 != UCS_OK)
    return status1;

  return UCS_OK;
}


ucs_status_t
ucx_mr_split_send_old(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
                  void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                  void *buffer1, ucs_memory_type_t mem_type1, int device1)
{
  // split_ratio : factor how much is sent to peer (0.39 -> 39% is sent to peer device)
  // element_size: nof bytes of splitted memory needs to be multiple of element_size
  // buf_size1   : maximum possible splittable memory

  cudaError_t c_status;
  ucs_status_t status;

  if (mem_type0 != UCS_MEMORY_TYPE_CUDA && mem_type1 != UCS_MEMORY_TYPE_CUDA)
  {
    printf("Split recv only implemented for CUDA Memory\n");
    return UCS_ERR_NOT_IMPLEMENTED;
  }

  size_t tmp;
  size_t msg_size1;
  size_t nof_elements;

  nof_elements = (size_t)roundf((msg_size0 / element_size) * (1 - split_ratio));

  DEBUG_PRINT("Full message size: %d Bytes, nof_elements: %d\n", msg_size0, nof_elements);

  tmp = nof_elements * element_size;
  msg_size1 = msg_size0 - tmp;
  msg_size0 = tmp;

  int canAccess = 0;

  DEBUG_PRINT("Calculated msg_size0: %d Bytes, msg_size_1: %d, Sum: %d\n", msg_size0, msg_size1, msg_size0 + msg_size1);

  DEBUG_PRINT("cudaMemcpyPeer: dst: %lld, dstDevice: %d, src: %lld, srcDevice: %d\n", buffer1, device1, buffer0 + msg_size0, device0);

  cudaDeviceCanAccessPeer(&canAccess, 0, 1);

  DEBUG_PRINT("Can Access: %d\n", canAccess);

  fflush(stdout);

  c_status = cudaMemcpyPeer(buffer1, device1, (buffer0 + msg_size0), device0, msg_size1);

  DEBUG_PRINT("CudaMemcpy Peer finished!\n");

  if (c_status != cudaSuccess)
    return UCS_ERR_REJECTED;

  status = ucx_mr_parallel_send(mr_ctx, tag,
                                buffer0, msg_size0, mem_type0, device0,
                                buffer1, msg_size1, mem_type1, device1);

  return status;
}



ucs_status_t
ucx_mr_single_send(ucx_mr_context_t *mr_ctx, int rail, ucp_tag_t tag,
                   void *buffer, size_t msg_size, ucs_memory_type_t mem_type, int device)
{
  ucs_status_t status;

  struct ucx_context *request = NULL;
  ucp_request_param_t send_param;

  ucp_ep_h *ep = mr_ctx->ep;
  ucp_worker_h *worker = mr_ctx->worker;

  ucx_mr_write_send_param(&send_param, mem_type);

  cudaSetDevice(rail);

  request = (struct ucx_context *)ucp_tag_send_nbx(ep[rail], buffer, msg_size, tag, &send_param);
  status = ucx_wait(worker[rail], request, "send", data_msg_str);

  return status;
}

ucs_status_t
ucx_mr_test_send(ucx_mr_context_t *mr_ctx, int rail, ucp_tag_t tag)
{
  ucs_status_t status;

  struct ucx_context *request = NULL;
  ucp_request_param_t send_param;

  ucp_ep_h ep = mr_ctx->ep[rail];
  ucp_worker_h worker = mr_ctx->worker[rail];

  size_t length = 10;

  void *buffer = malloc(length * sizeof(int));
  ucs_memory_type_t mem_type = UCS_MEMORY_TYPE_HOST;

  printf("Start to send!\n");
  ucx_mr_create_test_message(buffer, length, length * rail, mem_type);
  printf("Message: ");
  ucx_mr_read_test_message(buffer, length, mem_type);

  ucx_mr_write_send_param(&send_param, mem_type);

  request = (struct ucx_context *)ucp_tag_send_nbx(ep, buffer, length * sizeof(int), tag, &send_param);
  status = ucx_wait(worker, request, "send", data_msg_str);

  return status;
}


// // quad rail routines


// ucs_status_t
// ucx_mr_split_q_send(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
//                   void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
//                   void *buffer1, ucs_memory_type_t mem_type1, int device1,
//                   void *buffer2, ucs_memory_type_t mem_type2, int device2,
//                   void *buffer3, ucs_memory_type_t mem_type3, int device3)
// {
//   // split_ratio : factor how much is sent to peer (0.39 -> 39% is sent to peer device)
//   // element_size: nof bytes of splitted memory needs to be multiple of element_size
//   // buf_size1   : maximum possible splittable memory

//   cudaError_t c_status;
//   ucs_status_t status0;
//   ucs_status_t status1;

//   struct ucx_context *request0 = NULL;
//   struct ucx_context *request1 = NULL;
//   ucp_request_param_t send_param0;
//   ucp_request_param_t send_param1;

//   ucp_ep_h *ep = mr_ctx->ep;
//   ucp_worker_h *worker = mr_ctx->worker;

//   ucx_mr_write_send_param(&send_param0, mem_type0);
//   ucx_mr_write_send_param(&send_param1, mem_type1);

//   if (mem_type0 != UCS_MEMORY_TYPE_CUDA && mem_type1 != UCS_MEMORY_TYPE_CUDA)
//   {
//     printf("Split recv only implemented for CUDA Memory\n");
//     return UCS_ERR_NOT_IMPLEMENTED;
//   }

//   size_t tmp;
//   size_t msg_size1;
//   size_t nof_elements;

//   nof_elements = (size_t)roundf((msg_size0 / element_size) * (1 - split_ratio));

//   DEBUG_PRINT("Full message size: %d Bytes, nof_elements: %d\n", msg_size0, nof_elements);

//   tmp = nof_elements * element_size;
//   msg_size1 = msg_size0 - tmp;
//   msg_size0 = tmp;  

//   int canAccess = 0;

//   DEBUG_PRINT("Calculated msg_size0: %d Bytes, msg_size_1: %d, Sum: %d\n", msg_size0, msg_size1, msg_size0 + msg_size1);

//   DEBUG_PRINT("cudaMemcpyPeer: dst: %lld, dstDevice: %d, src: %lld, srcDevice: %d\n", buffer1, device1, buffer0 + msg_size0, device0);


//   cudaSetDevice(0);
//   request0 = (struct ucx_context *)ucp_tag_send_nbx(ep[0], buffer0, msg_size0, tag, &send_param0);
//    printf("s0\n");

//   c_status = cudaMemcpyPeer(buffer1, device1, (buffer0 + msg_size0), device0, msg_size1);
//   if (c_status != cudaSuccess)
//     return UCS_ERR_REJECTED;

//   DEBUG_PRINT("CudaMemcpy Peer finished!\n");

//   cudaSetDevice(1);
//   request1 = (struct ucx_context *)ucp_tag_send_nbx(ep[1], buffer1, msg_size1, tag, &send_param1);
//   printf("s1\n");




//   status0 = ucx_wait(worker[0], request0, "send", data_msg_str);
//   status1 = ucx_wait(worker[1], request1, "send", data_msg_str);

//   if (status0 != UCS_OK)
//     return status0;
//   if (status1 != UCS_OK)
//     return status1;

//   return UCS_OK;
// }




#endif // UCX_MR_SEND_H