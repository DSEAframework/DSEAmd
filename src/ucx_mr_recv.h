#ifndef UCX_MR_RECV_H
#define UCX_MR_RECV_H

#include <cuda_runtime.h>
#include <cuda.h>
#include <driver_types.h>

#include <math.h>

#include <ucp/api/ucp.h>

#include "ucx_mr_types.h"
#include "ucx_mr_wait.h"

static void
recv_handler(void *request, ucs_status_t status,
             const ucp_tag_recv_info_t *info, void *user_data)
{
  struct ucx_context *context = (struct ucx_context *)request;

  context->completed = 1;
}

void ucx_mr_write_recv_param(ucp_request_param_t *recv_param, ucs_memory_type_t mem_type)
{
  recv_param->op_attr_mask = UCP_OP_ATTR_FIELD_CALLBACK |
                             UCP_OP_ATTR_FIELD_DATATYPE |
                             UCP_OP_ATTR_FIELD_MEMORY_TYPE |
                             UCP_OP_ATTR_FLAG_NO_IMM_CMPL;

  recv_param->datatype = ucp_dt_make_contig(1);
  recv_param->memory_type = mem_type;
  recv_param->cb.recv = recv_handler;
}

ucs_status_t
ucx_mr_parallel_recv(ucx_mr_context_t *mr_ctx, ucp_tag_t tag,
                     void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                     void *buffer1, size_t msg_size1, ucs_memory_type_t mem_type1, int device1,
                     int * rec_size)
{
  ucs_status_t status0;
  ucs_status_t status1;

  struct ucx_context *request0 = NULL;
  struct ucx_context *request1 = NULL;
  ucp_tag_message_h msg_tag0 = NULL;
  ucp_tag_message_h msg_tag1 = NULL;
  ucp_tag_recv_info_t info_tag0;
  ucp_tag_recv_info_t info_tag1;

  ucp_request_param_t recv_param0;
  ucp_request_param_t recv_param1;

  ucp_ep_h *ep = mr_ctx->ep;
  ucp_worker_h *worker = mr_ctx->worker;

  ucx_mr_write_recv_param(&recv_param0, mem_type0);
  ucx_mr_write_recv_param(&recv_param1, mem_type1);

  for (;;)
  {
    if (msg_tag0 == NULL)
    {
      cudaSetDevice(device0);
      msg_tag0 = ucp_tag_probe_nb(worker[0], tag, tag_mask, 1, &info_tag0);

      if (msg_tag0 != NULL)
      {
        request0 = (struct ucx_context*)ucp_tag_msg_recv_nbx(worker[0], buffer0, msg_size0, msg_tag0, &recv_param0);
      }
      ucp_worker_progress(worker[0]);
    }

    if (msg_tag1 == NULL)
    {
      cudaSetDevice(device1);
      msg_tag1 = ucp_tag_probe_nb(worker[1], tag, tag_mask, 1, &info_tag1);

      if (msg_tag1 != NULL)
      {
        request1 = (ucx_context*)ucp_tag_msg_recv_nbx(worker[1], buffer1, msg_size1, msg_tag1, &recv_param1);
      }
      ucp_worker_progress(worker[1]);
    }

    if (msg_tag0 != NULL && msg_tag1 != NULL)
      break;
  }

  rec_size[0]=info_tag0.length;
  rec_size[1]=info_tag1.length;

  cudaSetDevice(device0);
  status0 = ucx_wait(worker[0], request0, "receive", data_msg_str);
  cudaSetDevice(device1);
  status1 = ucx_wait(worker[1], request1, "receive", data_msg_str);

  if (status0 != UCS_OK)
    return status0;
  if (status1 != UCS_OK)
    return status1;

  return UCS_OK;
}

ucs_status_t
ucx_mr_split_recv(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
                  void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                  void *buffer1, ucs_memory_type_t mem_type1, int device1)
{
  // 0 is the collecting device

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

  // tmp = nof_elements * element_size;
  // msg_size1 = msg_size0 - tmp;
  // msg_size0 = tmp;

  DEBUG_PRINT("Calculated msg_size0: %d Bytes, msg_size_1: %d, Sum: %d\n", msg_size0, msg_size1, msg_size0 + msg_size1);

  int rec_size[2]={};

  status = ucx_mr_parallel_recv(mr_ctx, tag,
                                buffer0, msg_size0, mem_type0, device0,
                                buffer1, msg_size0, mem_type1, device1,rec_size);

  if (status != UCS_OK)
    return status;

  // c_status = cudaMemcpyPeer((buffer0 + msg_size0), device0, buffer1, device1, msg_size1);
  c_status = cudaMemcpyPeer((buffer0 + rec_size[0]), device0, buffer1, device1, rec_size[1]);

  if (c_status != cudaSuccess)
    return UCS_ERR_REJECTED;

  return UCS_OK;
}

ucs_status_t
ucx_mr_single_recv(ucx_mr_context_t *mr_ctx, int rail, ucp_tag_t tag,
                   void *buffer, size_t msg_size, ucs_memory_type_t mem_type, int device)
{
  ucs_status_t status;

  struct ucx_context *request = NULL;
  ucp_tag_message_h msg_tag = NULL;
  ucp_tag_recv_info_t info_tag;

  ucp_request_param_t recv_param;

  ucp_ep_h *ep = mr_ctx->ep;
  ucp_worker_h *worker = mr_ctx->worker;

  ucx_mr_write_recv_param(&recv_param, mem_type);

  cudaSetDevice(rail);

  for (;;)
  {
    msg_tag = ucp_tag_probe_nb(worker[rail], tag, tag_mask, 1, &info_tag);
    if (msg_tag != NULL)
    {
      break;
    }

    ucp_worker_progress(worker[rail]);
  }

  request = (struct ucx_context *)ucp_tag_msg_recv_nbx(worker[rail], buffer, msg_size, msg_tag, &recv_param);

  status = ucx_wait(worker[rail], request, "receive", data_msg_str);

  return status;
}

ucs_status_t
ucx_mr_test_recv(ucx_mr_context_t *mr_ctx, int rail, ucp_tag_t tag)
{
  ucs_status_t status;

  struct ucx_context *request = NULL;
  ucp_tag_message_h msg_tag;
  ucp_tag_recv_info_t info_tag;

  ucp_request_param_t recv_param;

  ucp_ep_h ep = mr_ctx->ep[rail];
  ucp_worker_h worker = mr_ctx->worker[rail];

  size_t length = 10;

  void *buffer = malloc(length * sizeof(int));
  ucs_memory_type_t mem_type = UCS_MEMORY_TYPE_HOST;

  ucx_mr_write_recv_param(&recv_param, mem_type);

  printf("Start to receive!\n");
  ucx_mr_init_test_message(buffer, length, 0, mem_type);
  printf("Message previous: ");
  ucx_mr_read_test_message(buffer, length, mem_type);

  for (;;)
  {
    msg_tag = ucp_tag_probe_nb(worker, tag, tag_mask, 1, &info_tag);

    if (msg_tag != NULL)
    {
      break;
    }
    ucp_worker_progress(worker);
  }

  printf("Message is here! \n");

  request = (struct ucx_context *)ucp_tag_msg_recv_nbx(worker, buffer, length * sizeof(int), msg_tag, &recv_param);

  status = ucx_wait(worker, request, "receive", data_msg_str);

  printf("Message received: ");
  ucx_mr_read_test_message(buffer, length, mem_type);

  return status;
}











// qual rail routines


ucs_status_t
ucx_mr_parallel_q_recv(ucx_mr_context_t *mr_ctx, ucp_tag_t tag,
                     void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                     void *buffer1, size_t msg_size1, ucs_memory_type_t mem_type1, int device1,
                     void *buffer2, size_t msg_size2, ucs_memory_type_t mem_type2, int device2,
                     void *buffer3, size_t msg_size3, ucs_memory_type_t mem_type3, int device3,
                     int * rec_size)
{
  ucs_status_t status0;
  ucs_status_t status1;
  ucs_status_t status2;
  ucs_status_t status3;

  struct ucx_context *request0 = NULL;
  struct ucx_context *request1 = NULL;
  struct ucx_context *request2 = NULL;
  struct ucx_context *request3 = NULL;

  ucp_tag_message_h msg_tag0 = NULL;
  ucp_tag_message_h msg_tag1 = NULL;
  ucp_tag_message_h msg_tag2 = NULL;
  ucp_tag_message_h msg_tag3 = NULL;
  ucp_tag_recv_info_t info_tag0;
  ucp_tag_recv_info_t info_tag1;
  ucp_tag_recv_info_t info_tag2;
  ucp_tag_recv_info_t info_tag3;

  ucp_request_param_t recv_param0;
  ucp_request_param_t recv_param1;
  ucp_request_param_t recv_param2;
  ucp_request_param_t recv_param3;

  ucp_ep_h *ep = mr_ctx->ep;
  ucp_worker_h *worker = mr_ctx->worker;

  ucx_mr_write_recv_param(&recv_param0, mem_type0);
  ucx_mr_write_recv_param(&recv_param1, mem_type1);
  ucx_mr_write_recv_param(&recv_param2, mem_type2);
  ucx_mr_write_recv_param(&recv_param3, mem_type3);

  for (;;)
  {
    if (msg_tag0 == NULL) {
      cudaSetDevice(device0);
      msg_tag0 = ucp_tag_probe_nb(worker[0], tag, tag_mask, 1, &info_tag0);

      if (msg_tag0 != NULL) {
        request0 = (struct ucx_context*)ucp_tag_msg_recv_nbx(worker[0], buffer0, msg_size0, msg_tag0, &recv_param0);
      }
      ucp_worker_progress(worker[0]);
    }

    if (msg_tag1 == NULL) {
      cudaSetDevice(device1);
      msg_tag1 = ucp_tag_probe_nb(worker[1], tag, tag_mask, 1, &info_tag1);

      if (msg_tag1 != NULL) {
        request1 = (ucx_context*)ucp_tag_msg_recv_nbx(worker[1], buffer1, msg_size1, msg_tag1, &recv_param1);
      }
      ucp_worker_progress(worker[1]);
    }

    if (msg_tag2 == NULL) {
      cudaSetDevice(device2);
      msg_tag2 = ucp_tag_probe_nb(worker[2], tag, tag_mask, 1, &info_tag2);

      if (msg_tag2 != NULL) {
        request2 = (ucx_context*)ucp_tag_msg_recv_nbx(worker[2], buffer2, msg_size2, msg_tag2, &recv_param2);
      }
      ucp_worker_progress(worker[2]);
    }

    if (msg_tag3 == NULL) {
      cudaSetDevice(device3);
      msg_tag3 = ucp_tag_probe_nb(worker[3], tag, tag_mask, 1, &info_tag3);

      if (msg_tag3 != NULL) {
        request3 = (ucx_context*)ucp_tag_msg_recv_nbx(worker[3], buffer3, msg_size3, msg_tag3, &recv_param3);
      }
      ucp_worker_progress(worker[3]);
    }

    if (msg_tag0 != NULL && msg_tag1 != NULL && msg_tag2 != NULL && msg_tag3 != NULL)
      break;
  }

  rec_size[0]=info_tag0.length;
  rec_size[1]=info_tag1.length;
  rec_size[2]=info_tag2.length;
  rec_size[3]=info_tag3.length;

  cudaSetDevice(device0);
  status0 = ucx_wait(worker[0], request0, "receive", data_msg_str);
  printf("r0\n");
  cudaSetDevice(device1);
  status1 = ucx_wait(worker[1], request1, "receive", data_msg_str);
  printf("r1\n");
  cudaSetDevice(device2);
  status1 = ucx_wait(worker[2], request2, "receive", data_msg_str);
  printf("r2\n");
  cudaSetDevice(device3);
  status1 = ucx_wait(worker[3], request3, "receive", data_msg_str);
  printf("r3\n");

  if (status0 != UCS_OK)
    return status0;
  if (status1 != UCS_OK)
    return status1;
  if (status2 != UCS_OK)
    return status2;
  if (status3 != UCS_OK)
    return status3;



  return UCS_OK;
}



ucs_status_t
ucx_mr_split_q_recv(ucx_mr_context_t *mr_ctx, ucp_tag_t tag, float split_ratio, int element_size,
                  void *buffer0, size_t msg_size0, ucs_memory_type_t mem_type0, int device0,
                  void *buffer1, ucs_memory_type_t mem_type1, int device1,
                  void *buffer2, ucs_memory_type_t mem_type2, int device2,
                  void *buffer3, ucs_memory_type_t mem_type3, int device3)
{
  // 0 is the collecting device

  cudaError_t c_status;
  ucs_status_t status;

  if (mem_type0 != UCS_MEMORY_TYPE_CUDA && mem_type1 != UCS_MEMORY_TYPE_CUDA)
  {
    printf("Split recv only implemented for CUDA Memory\n");
    return UCS_ERR_NOT_IMPLEMENTED;
  }
  if (mem_type2 != UCS_MEMORY_TYPE_CUDA && mem_type3 != UCS_MEMORY_TYPE_CUDA)
  {
    printf("Split recv only implemented for CUDA Memory\n");
    return UCS_ERR_NOT_IMPLEMENTED;
  }

  size_t tmp;
  size_t msg_size1;
  size_t nof_elements;

  nof_elements = (size_t)roundf((msg_size0 / element_size) * (1 - split_ratio));

  DEBUG_PRINT("Full message size: %d Bytes, nof_elements: %d\n", msg_size0, nof_elements);

  // tmp = nof_elements * element_size;
  // msg_size1 = msg_size0 - tmp;
  // msg_size0 = tmp;

  DEBUG_PRINT("Calculated msg_size0: %d Bytes, msg_size_1: %d, Sum: %d\n", msg_size0, msg_size1, msg_size0 + msg_size1);

  int rec_size[4]={};

  status = ucx_mr_parallel_q_recv(mr_ctx, tag,
                                buffer0, msg_size0, mem_type0, device0,
                                buffer1, msg_size0, mem_type1, device1,
                                buffer2, msg_size0, mem_type2, device2,
                                buffer3, msg_size0, mem_type3, device3,
                                rec_size);

  if (status != UCS_OK)
    return status;

  // c_status = cudaMemcpyPeer((buffer0 + msg_size0), device0, buffer1, device1, msg_size1);
  c_status = cudaMemcpyPeer((buffer0 + rec_size[0]), device0, buffer1, device1, rec_size[1]);
// todo!!!


  if (c_status != cudaSuccess)
    return UCS_ERR_REJECTED;

  return UCS_OK;
}




#endif // UCX_MR_RECV_H