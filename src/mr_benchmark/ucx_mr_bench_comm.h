#ifndef UCX_MR_BENCH_COMM_H
#define UCX_MR_BENCH_COMM_H

#include <sys/times.h>
#include <unistd.h>

#include <omp.h>

#include <ucp/api/ucp.h>

#include "ucx_mr_bench_types.h"
#include "../ucx_multirail.h"

#ifndef MAX_RUNS
#define MAX_RUNS 500
#endif


void ucx_mr_bench_test_send_split(ucx_mr_bench_context_t *mr_bench_ctx, float split_ratio)
{
  ucs_status_t status;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;

  void** buffer = mr_bench_ctx->send_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;


  int element_size = 4;


  printf("Start to send!\n");
  ucx_mr_create_test_message(buffer[0], length/element_size, 1, mem_type);
  printf("Message: ");
  ucx_mr_read_test_message(buffer[0], length/element_size, mem_type);

  

  ucx_mr_split_send(mr_ctx, tag, split_ratio, element_size,
                      buffer[0], length, mem_type, 0,
                      buffer[1], mem_type, 1);

}

void ucx_mr_bench_test_recv_split(ucx_mr_bench_context_t *mr_bench_ctx, float split_ratio)
{
  ucs_status_t status;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;

  void** buffer = mr_bench_ctx->recv_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;


  int element_size = 4;

  printf("Before Reiving: ");
  ucx_mr_read_test_message(buffer[0], length/element_size, mem_type);


  ucx_mr_split_recv(mr_ctx, tag, split_ratio, element_size,
                      buffer[0], length, mem_type, 0,
                      buffer[1], mem_type, 1);


  printf("Received Message: ");
  ucx_mr_read_test_message(buffer[0], length/element_size, mem_type);

}


void ucx_mr_bench_send_split_old(ucx_mr_bench_context_t *mr_bench_ctx, float split_ratio)
{
  ucs_status_t status;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void **buffer = mr_bench_ctx->send_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int element_size = 4;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < 2 * MAX_RUNS; ++runs)
  {

    ucx_mr_split_send_old(mr_ctx, tag, split_ratio, element_size,
                      buffer[0], length, mem_type, 0,
                      buffer[1], mem_type, 1);

    DEBUG_PRINT("Run Finished!\n\n\n");
  }
  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}


void ucx_mr_bench_send_split(ucx_mr_bench_context_t *mr_bench_ctx, float split_ratio)
{
  ucs_status_t status;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void **buffer = mr_bench_ctx->send_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int element_size = 4;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < 2 * MAX_RUNS; ++runs)
  {

    ucx_mr_split_send(mr_ctx, tag, split_ratio, element_size,
                      buffer[0], length, mem_type, 0,
                      buffer[1], mem_type, 1);

    DEBUG_PRINT("Run Finished!\n\n\n");
  }
  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}

void ucx_mr_bench_recv_split(ucx_mr_bench_context_t *mr_bench_ctx, float split_ratio)
{
  ucs_status_t status;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void **buffer = mr_bench_ctx->recv_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int element_size = 4;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < 2 * MAX_RUNS; ++runs)
  {
    DEBUG_PRINT("Run: %d\n", runs);

    ucx_mr_split_recv(mr_ctx, tag, split_ratio, element_size,
                      buffer[0], length, mem_type, 0,
                      buffer[1], mem_type, 1);

    DEBUG_PRINT("Run Finished!\n\n\n");
  }

  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}

void ucx_mr_bench_send_dual_mp(ucx_mr_bench_context_t *mr_bench_ctx)
{
  ucs_status_t status0;
  ucs_status_t status1;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void **buffer = mr_bench_ctx->send_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < MAX_RUNS; ++runs)
  {
    DEBUG_PRINT("Run: %d\n", runs);
#pragma omp parallel
    {
#pragma omp single nowait
      {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        DEBUG_PRINT("Send message with tag: %llx, from thread: %d\n", tag, t);

        ucx_mr_single_send(mr_ctx, 0, tag, buffer[0], length, mem_type, 0);
      }
#pragma omp single nowait
      {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        DEBUG_PRINT("Send message with tag: %llx, from thread: %d\n", tag, t);

        ucx_mr_single_send(mr_ctx, 1, tag, buffer[1], length, mem_type, 10);
      }
    }
    /* pragma parallel end */
    DEBUG_PRINT("Run Finished!\n\n\n");
  }

  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: 2 * %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}

void ucx_mr_bench_recv_dual_mp(ucx_mr_bench_context_t *mr_bench_ctx)
{
  ucs_status_t status;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void **buffer = mr_bench_ctx->send_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < MAX_RUNS; ++runs)
  {
    DEBUG_PRINT("Run: %d\n", runs);

#pragma omp parallel
    {
#pragma omp single nowait
      {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        ucx_mr_single_recv(mr_ctx, 0, tag, buffer[0], length, mem_type, 0);
        DEBUG_PRINT("End Task 0\n");
      }
#pragma omp single nowait
      {
        int t = 0;
#ifdef _OPENMP
        t = omp_get_thread_num();
#endif
        ucx_mr_single_recv(mr_ctx, 1, tag, buffer[1], length, mem_type, 1);
        DEBUG_PRINT("End Task 1\n");
      }
    }
    /* omp end parallel*/
    DEBUG_PRINT("Run Finished!\n\n\n");
  }

  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: 2 * %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}

void ucx_mr_bench_send_dual(ucx_mr_bench_context_t *mr_bench_ctx)
{
  ucs_status_t status;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void **buffer = mr_bench_ctx->send_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < MAX_RUNS; ++runs)
  {

    ucx_mr_parallel_send(mr_ctx, tag,
                         buffer[0], length, mem_type, 0,
                         buffer[1], length, mem_type, 1);

    DEBUG_PRINT("Run Finished!\n\n\n");
  }
  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: 2 * %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}

void ucx_mr_bench_recv_dual(ucx_mr_bench_context_t *mr_bench_ctx)
{
  ucs_status_t status;

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void **buffer = mr_bench_ctx->recv_buffer;
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < MAX_RUNS; ++runs)
  {
    DEBUG_PRINT("Run: %d\n", runs);

    status = ucx_mr_parallel_recv(mr_ctx, tag,
                                  buffer[0], length, mem_type, 0,
                                  buffer[1], length, mem_type, 1);

    DEBUG_PRINT("Run Finished!\n\n\n");
  }

  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: 2 * %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}

void ucx_mr_bench_recv_single(ucx_mr_bench_context_t *mr_bench_ctx, int rail)
{
  ucs_status_t status;

  cudaSetDevice(rail);

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void *buffer = mr_bench_ctx->recv_buffer[rail];
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < MAX_RUNS; ++runs)
  {

    DEBUG_PRINT("Run: %d\n", runs);

    status = ucx_mr_single_recv(mr_ctx, rail, tag, buffer, length, mem_type, rail);
  }

  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}

void ucx_mr_bench_send_single(ucx_mr_bench_context_t *mr_bench_ctx, int rail)
{
  ucs_status_t status;

  cudaSetDevice(rail);

  ucx_mr_context_t *mr_ctx = &mr_bench_ctx->mr_ctx;
  ucp_tag_t tag = 0x33;
  void *buffer = mr_bench_ctx->send_buffer[rail];
  ucs_memory_type_t mem_type = mr_bench_ctx->mem_type;
  size_t length = mr_bench_ctx->msg_size;

  int runs;
  double t0 = wallTime();
  for (runs = 0; runs < MAX_RUNS; ++runs)
  {
    DEBUG_PRINT("Run: %d Send message with tag: %llx\n", runs, tag + runs);

    ucx_mr_single_send(mr_ctx, rail, tag, buffer, length, mem_type, rail);

    DEBUG_PRINT("Run Finished!\n\n\n");
  }
  double t1 = wallTime() - t0;
  printf("\n\n\n====================\n");

  printf("%d runs took %lfs to send.\n", runs, t1);

  double msg_size = length * 1e-6;
  double t_per_run = t1 / runs;
  printf("Msg size: %lf MB, Bandwidth: %lf MB/s\n", msg_size, msg_size / t_per_run);

  printf("====================\n\n\n");
}

#endif // UCX_MR_BENCH_COMM_H