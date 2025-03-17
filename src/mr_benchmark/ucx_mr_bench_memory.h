#ifndef MR_BENCHMARK_UCX_MR_BENCH_MEMORY_H
#define MR_BENCHMARK_UCX_MR_BENCH_MEMORY_H

void ucx_mr_free(void *buffer, ucs_memory_type_t mem_type)
{
  if (mem_type == UCS_MEMORY_TYPE_HOST)
  {
    free(buffer);
  }
  else if (mem_type == UCS_MEMORY_TYPE_CUDA)
  {
    cudaFree(buffer);
  }
}

void ucx_mr_free_mem(ucx_mr_bench_context_t *mr_ctx)
{
  DEBUG_PRINT("Free Memory.\n");
  for (int i = 0; i < mr_ctx->recv_buffer_count; ++i)
  {
    ucx_mr_free(mr_ctx->recv_buffer[i], mr_ctx->mem_type);
  }
  for (int i = 0; i < mr_ctx->send_buffer_count; ++i)
  {
    ucx_mr_free(mr_ctx->send_buffer[i], mr_ctx->mem_type);
  }
}

ucs_status_t
ucx_mr_alloc(void **buffer, size_t length, ucs_memory_type_t mem_type, int device)
{
  ucs_status_t status = UCS_OK;
  if (mem_type == UCS_MEMORY_TYPE_HOST)
  {
    *buffer = malloc(length);
    if (*buffer == NULL)
    {
      ucs_error("failed to allocate host memory");
      status = UCS_ERR_NO_MEMORY;
    }
  }
  else if (mem_type == UCS_MEMORY_TYPE_CUDA)
  {
    cudaSetDevice(device);
    cudaError_t e = cudaMalloc(buffer, length);
    if (e != cudaSuccess)
    {
      printf("Error in CudaMalloc.\n\n");
      status = UCS_ERR_NO_MEMORY;
    }
  }

  return status;
}

ucs_status_t
ucx_mr_alloc_mem(ucx_mr_bench_context_t *mr_ctx)
{
  ucs_status_t status;

  for (int i = 0; i < NOF_RAILS; ++i)
  {
    int device = i % NOF_RAILS;
    status = ucx_mr_alloc(&mr_ctx->recv_buffer[i], mr_ctx->msg_size, mr_ctx->mem_type, device);
    if (status != UCS_OK)
    {
      goto err_free_mem;
    }
    mr_ctx->recv_buffer_count = i + 1;
  }

  for (int i = 0; i < NOF_RAILS; ++i)
  {
    int device = i % NOF_RAILS;
    status = ucx_mr_alloc(&mr_ctx->send_buffer[i], mr_ctx->msg_size, mr_ctx->mem_type, device);
    if (status != UCS_OK)
    {
      goto err_free_mem;
    }
    mr_ctx->send_buffer_count = i + 1;
  }

  return status;

err_free_mem:
  ucx_mr_free_mem(mr_ctx);
err:
  return status;
}

#endif // MR_BENCHMARK_UCX_MR_BENCH_MEMORY_H
