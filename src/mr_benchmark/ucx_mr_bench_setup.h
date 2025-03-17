#ifndef MR_BENCHMARK_UCX_MR_BENCH_SETUP_H
#define MR_BENCHMARK_UCX_MR_BENCH_SETUP_H

#include <ucp/api/ucp.h>
#include "ucx_mr_bench_types.h"
#include "ucx_mr_bench_memory.h"

ucs_status_t
ucx_mr_bench_setup(ucx_mr_bench_context_t *mr_bench_ctx)
{
  ucs_status_t status;
  mr_bench_ctx->recv_buffer_count = 0;
  mr_bench_ctx->send_buffer_count = 0;
  mr_bench_ctx->repeated_test = 0;

  printf("Setup UCX!\n");
  status = ucx_mr_setup(&mr_bench_ctx->mr_ctx);
  if (status != UCS_OK)
  {
    printf("Something went wrong!\n");
    return status;
  }

  printf("Allocate Memory !\n");
  status = ucx_mr_alloc_mem(mr_bench_ctx);

  return status;
}

#endif // MR_BENCHMARK_UCX_MR_BENCH_SETUP_H
