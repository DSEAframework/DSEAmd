#ifndef MR_BENCHMARK_UCX_MR_BENCH_H
#define MR_BENCHMARK_UCX_MR_BENCH_H

#include <ucp/api/ucp.h>
#include <unistd.h>
#include <string.h>
#include "ucx_mr_bench_types.h"
#include "ucx_mr_bench_setup.h"
#include "ucx_mr_bench_memory.h"
#include "ucx_mr_bench_comm.h"

#include "../ucx_mr_aux.h"

ucs_status_t
parse_bench_opts(ucx_mr_bench_context_t *mr_bench_ctx, int argc, char **argv)
{
    ucs_status_t status;
    int c;
    char *ptr;

    ucs_trace_func("");

    optind = 1;
    while ((c = getopt(argc, argv, "p:R:A:T:M:s:P:")) != -1)
    {
        switch (c)
        {
        case 'T':
            DEBUG_PRINT("Got Test %s\n", optarg);
            if (!strcmp(optarg, "DUAL"))
            {
                mr_bench_ctx->test_type = DUAL;
            }
            else if (!strcmp(optarg, "DUAL_MP"))
            {
                mr_bench_ctx->test_type = DUAL_MP;
            }
            else if (!strcmp(optarg, "SINGLE"))
            {
                mr_bench_ctx->test_type = SINGLE;
            }
            else if (!strcmp(optarg, "SPLIT"))
            {
                mr_bench_ctx->test_type = SPLIT;
            }
            else if (!strcmp(optarg, "TEST_SPLIT"))
            {
                mr_bench_ctx->test_type = TEST_SPLIT;
            }
            else
            {
                return UCS_ERR_INVALID_PARAM;
            }
            break;
        case 'M':
            DEBUG_PRINT("Got Memory %s\n", optarg);
            if (!strcmp(optarg, "CUDA"))
            {
                mr_bench_ctx->mem_type = UCS_MEMORY_TYPE_CUDA;
            }
            else if (!strcmp(optarg, "HOST"))
            {
                mr_bench_ctx->mem_type = UCS_MEMORY_TYPE_HOST;
            }
            else
            {
                return UCS_ERR_INVALID_PARAM;
            }
            break;
        case 's':
            DEBUG_PRINT("Got Message size %s\n", optarg);
            mr_bench_ctx->msg_size = atoi(optarg);
            break;
        case 'p':
            DEBUG_PRINT("Got a port %s\n", optarg);
            mr_bench_ctx->mr_ctx.port = atoi(optarg);
            break;
        case 'R':
            ptr = strtok(optarg, ",");
            mr_bench_ctx->mr_ctx.rail0 = ptr;
            DEBUG_PRINT("Got Rail %s\n", ptr);
            ptr = strtok(NULL, ",");
            mr_bench_ctx->mr_ctx.rail1 = ptr;
            DEBUG_PRINT("Got Rail %s\n", ptr);
            break;
        case 'A':
            DEBUG_PRINT("Got Server Address %s\n", optarg);
            mr_bench_ctx->mr_ctx.server_addr = optarg;
            break;
        case 'P':
            DEBUG_PRINT("Got Repeated test with %s repetitions per msg size\n", optarg);
            mr_bench_ctx->repetitions = atoi(optarg);
            mr_bench_ctx->repeated_test = 1;
            break;
        default:
            DEBUG_PRINT("Default\n");
            break;
        }
    }

    return UCS_OK;
}

#endif // MR_BENCHMARK_UCX_MR_BENCH_H