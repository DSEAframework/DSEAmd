#ifndef UCX_MR_CLEANUP
#define UCX_MR_CLEANUP

#include <ucp/api/ucp.h>

#include "ucx_mr_types.h"
#include "ucx_mr_sock_comm.h"

void ucx_mr_destroy_endpoint(ucx_mr_context_t *mr_ctx)
{
  for (int i = 0; i < mr_ctx->ep_count; ++i)
  {
    // Not ready
    ucs_status_ptr_t req;
    ucs_status_t status;

    ucx_mr_worker_progress(mr_ctx, i);

    if (mr_ctx->ep[i] != NULL)
    {
      req = ucp_ep_close_nb(mr_ctx->ep[i], UCP_EP_CLOSE_MODE_FLUSH);
    }
    if (!UCS_PTR_IS_PTR(req) && req != NULL)
    {
      ucs_warn("failed to close ep %p: %s\n",
               mr_ctx->ep[i],
               ucs_status_string(UCS_PTR_STATUS(req)));
    }

    ucx_mr_worker_progress(mr_ctx, i);
    if (req)
    {
      status = ucp_request_check_status(req);
      while (status == UCS_INPROGRESS)
      {
        ucx_mr_worker_progress(mr_ctx, i);
        status = ucp_request_check_status(req);
      }
      ucp_request_free(req);
    }
  }
}

void ucx_mr_destroy_worker(ucx_mr_context_t *mr_ctx)
{
  for (int i = 0; i < mr_ctx->worker_count; ++i)
  {
    ucp_worker_destroy(mr_ctx->worker[i]);
  }
}

void ucx_mr_destroy_ctx(ucx_mr_context_t *mr_ctx)
{
  for (int i = 0; i < mr_ctx->ctx_count; ++i)
  {
    ucp_cleanup(mr_ctx->ctx[i]);
  };
}

void ucx_mr_cleanup(ucx_mr_context_t *mr_ctx, CleanUp pos)
{
  switch (pos)
  {
  case FULL:
  case EP:
    ucx_mr_destroy_endpoint(mr_ctx);
    DEBUG_PRINT("Endpoints destroyed.\n");
  case WORKER:
    ucx_mr_destroy_worker(mr_ctx);
    DEBUG_PRINT("Worker destroyed.\n");
  case CTX:
    ucx_mr_destroy_ctx(mr_ctx);
    DEBUG_PRINT("UCP Cleanup.\n");
  case COMM:
    ucx_mr_cleanup_comm(mr_ctx);
    DEBUG_PRINT("Comm cleanup.\n");
  default:;
  }
}

#endif // UCX_MR_CLEANUP