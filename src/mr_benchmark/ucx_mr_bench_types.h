#ifndef MR_BENCHMARK_UCX_MR_BENCH_TYPES_H
#define MR_BENCHMARK_UCX_MR_BENCH_TYPES_H

#include "../ucx_multirail.h"

typedef enum MR_Benchmarks
{
    DUAL,
    DUAL_MP,
    SINGLE,
    SPLIT,
    TEST_SPLIT
} MR_Benchmarks;

typedef struct ucx_mr_bench_context
{

    // Multirail context
    ucx_mr_context_t mr_ctx;
    // Memory
    size_t msg_size;
    int recv_buffer_count;
    int send_buffer_count;
    void *recv_buffer[NOF_RAILS];
    void *send_buffer[NOF_RAILS];

    ucs_memory_type_t mem_type;

    MR_Benchmarks test_type;
    int           repetitions;
    int           repeated_test;

} ucx_mr_bench_context_t;

#endif // MR_BENCHMARK_UCX_MR_BENCH_TYPES_H