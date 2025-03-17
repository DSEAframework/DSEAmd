#ifndef UCX_MR_TYPES_H
#define UCX_MR_TYPES_H

#include <ucp/api/ucp.h>

#define ADDR_BUF_SIZE 4096

#define NOF_RAILS 4

static const ucp_tag_t tag_mask = UINT64_MAX;
static const char *data_msg_str = "UCX data message";

typedef struct
{
    int num_outstanding; /* Number of outstanding flush operations */
    ucs_status_t status; /* Cumulative status of all flush operations */
} ucp_mr_flush_context_t;

struct ucx_context
{
    int completed;
};

typedef enum CommRole
{
    SENDER,
    RECIEVER
} CommRole;

typedef enum CleanUp
{
    FULL,
    CTX,
    COMM,
    WORKER,
    EP
} CleanUp;

typedef struct sock_rte_group
{
    int sendfd;
    int recvfd;
    int is_server;
    int size;
    int peer;
} sock_rte_group_t;

typedef struct ucx_mr_context
{
    // ucx objects
    ucp_context_h ctx[NOF_RAILS];
    ucp_worker_h worker[NOF_RAILS];
    ucp_ep_h ep[NOF_RAILS];

    // Auxiliary
    int ctx_count;
    int worker_count;
    int ep_count;
    CommRole role;

    const char *rail0;
    const char *rail1;
    const char *rail2;
    const char *rail3;

    // socket communictation
    const char *server_addr;
    sock_rte_group_t sock_rte_group;
    int port;

} ucx_mr_context_t;

#endif // UCX_MR_TYPES_H
