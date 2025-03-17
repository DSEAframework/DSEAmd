#ifndef UCX_MR_AUX_H
#define UCX_MR_AUX_H

#include <cuda_runtime.h>
#include <cuda.h>

#include <ucp/api/ucp.h>
#include <ucs/sys/sys.h>
#include <ucs/debug/log.h>

#include <sys/times.h>
#include <unistd.h>

#include "ucx_mr_types.h"

// #define DEBUG

#ifdef DEBUG
#define DEBUG_PRINT(fmt, args...) fprintf(stdout, "DEBUG: " fmt, ##args)
#else
#define DEBUG_PRINT(fmt, args...) /* Don't do anything in release builds */
#endif

double wallTime()
{
    static int ticks_per_second = 0;
    if (!ticks_per_second)
    {
        ticks_per_second = sysconf(_SC_CLK_TCK);
    }
    struct tms timebuf;
    /* times returns the number of real time ticks passed since start */
    return (double)times(&timebuf) / ticks_per_second;
}

static void
ucx_mr_init_test_message(void *buffer, size_t write_length, int init_val, ucs_memory_type_t mem_type)
{

    int arr[write_length];
    for (int i = 0; i < write_length; ++i)
    {
        arr[i] = init_val;
    }

    if (mem_type == UCS_MEMORY_TYPE_CUDA)
    {
        cudaMemcpy(buffer, arr, write_length * sizeof(int), cudaMemcpyHostToDevice);
    }
    else if (mem_type == UCS_MEMORY_TYPE_HOST)
    {
        for (int i = 0; i < write_length; ++i)
        {
            ((int *)buffer)[i] = arr[i];
        }
    }
}

static void
ucx_mr_create_test_message(void *buffer, size_t write_length, int start_val, ucs_memory_type_t mem_type)
{

    int arr[write_length];
    for (int i = 0; i < write_length; ++i)
    {
        arr[i] = start_val + i;
    }

    if (mem_type == UCS_MEMORY_TYPE_CUDA)
    {
        cudaMemcpy(buffer, arr, write_length * sizeof(int), cudaMemcpyHostToDevice);
    }
    else if (mem_type == UCS_MEMORY_TYPE_HOST)
    {
        for (int i = 0; i < write_length; ++i)
        {
            ((int *)buffer)[i] = arr[i];
        }
    }
}

static void
ucx_mr_read_test_message(void *buffer, size_t read_length, ucs_memory_type_t mem_type)
{
    int arr[read_length];

    if (mem_type == UCS_MEMORY_TYPE_CUDA)
    {
        cudaMemcpy(arr, buffer, read_length * sizeof(int), cudaMemcpyDeviceToHost);
    }
    else if (mem_type == UCS_MEMORY_TYPE_HOST)
    {
        for (int i = 0; i < read_length; ++i)
        {
            arr[i] = ((int *)buffer)[i];
        }
    }

    for (int i = 0; i < read_length; ++i)
    {
        printf("%d, ", arr[i]);
    }
    printf("\n");
}

static void
ucp_mr_worker_flush_callback(void *request, ucs_status_t status, void *user_data)
{
    ucp_mr_flush_context_t *ctx = (ucp_mr_flush_context_t*)user_data;

    --ctx->num_outstanding;
    if (status != UCS_OK)
    {
        ucs_error("worker flush callback got status %s",
                  ucs_status_string(status));
        ctx->status = status;
    }
    ucp_request_free(request);
}

static ucs_status_t
ucx_mr_flush_workers(ucx_mr_context_t *mr_ctx)
{
    ucp_mr_flush_context_t ctx = {
        .num_outstanding = 0,
        .status = UCS_OK};
    // ucp_request_param_t param = {
    //     .op_attr_mask = UCP_OP_ATTR_FIELD_CALLBACK |
    //                     UCP_OP_ATTR_FIELD_USER_DATA,
    //     .cb.send = ucp_mr_worker_flush_callback,
    //     .user_data = &ctx};

    ucp_request_param_t param;
    
        param.op_attr_mask = UCP_OP_ATTR_FIELD_CALLBACK |UCP_OP_ATTR_FIELD_USER_DATA;
        param.cb.send = ucp_mr_worker_flush_callback;
        param.user_data = &ctx;

    void *flush_req;
    unsigned i;

    for (i = 0; i < mr_ctx->worker_count; ++i)
    {
        flush_req = ucp_worker_flush_nbx(mr_ctx->worker[i], &param);
        if (UCS_PTR_IS_ERR(flush_req))
        {
            ctx.status = UCS_PTR_STATUS(flush_req);
            ucs_error("ucp_worker_flush_nbx() failed on thread %d: %s", i,
                      ucs_status_string(ctx.status));
        }

        if (UCS_PTR_IS_PTR(flush_req))
        {
            ++ctx.num_outstanding;
        }
    }

    /* Progress all workers in parallel to avoid deadlocks */
    while (ctx.num_outstanding > 0)
    {
        for (i = 0; i < mr_ctx->worker_count; ++i)
        {
            ucp_worker_progress(mr_ctx->worker[i]);
        }
    }

    return ctx.status;
}

static void
ucx_mr_worker_progress(void *arg, int idx)
{
    ucx_mr_context_t *mr_ctx = (ucx_mr_context_t *)arg;
    ucp_worker_progress(mr_ctx->worker[idx]);
}

static void
request_init(void *request)
{
    struct ucx_context *contex = (struct ucx_context *)request;

    contex->completed = 0;
}

ucs_status_t
parse_opts(ucx_mr_context_t *mr_ctx, int argc, char **argv)
{
    ucs_status_t status;
    int c;
    char *ptr;

    ucs_trace_func("");

    optind = 1;
    while ((c = getopt(argc, argv, "p:R:A:")) != -1)
    {
        switch (c)
        {
        case 'p':
            DEBUG_PRINT("Got a port %s\n", optarg);
            mr_ctx->port = atoi(optarg);
            break;
        case 'R':
            ptr = strtok(optarg, ",");
            mr_ctx->rail0 = ptr;
            DEBUG_PRINT("Got Rail %s\n", ptr);
            ptr = strtok(NULL, ",");
            mr_ctx->rail1 = ptr;
            DEBUG_PRINT("Got Rail %s\n", ptr);
            ptr = strtok(NULL, ",");
            mr_ctx->rail2 = ptr;
            DEBUG_PRINT("Got Rail %s\n", ptr);
            ptr = strtok(NULL, ",");
            mr_ctx->rail3 = ptr;
            DEBUG_PRINT("Got Rail %s\n", ptr);
            break;
        case 'A':
            DEBUG_PRINT("Got Server Address %s\n", optarg);
            mr_ctx->server_addr = optarg;
            break;
        default:
            DEBUG_PRINT("Default\n");
            break;
        }
    }

    return UCS_OK;
}

#endif // UCX_MR_AUX_H
