#ifndef UCX_MR_SETUP_H
#define UCX_MR_SETUP_H

#include <cuda_runtime.h>
#include <cuda.h>

#include <ucp/api/ucp.h>

#include "ucx_mr_types.h"
#include "ucx_mr_aux.h"
#include "ucx_mr_sock_comm.h"
#include "ucx_mr_cleanup.h"

#include "ucx_mr_recv.h"
#include "ucx_mr_send.h"

void ucx_mr_test_connection(ucx_mr_context_t *mr_ctx)
{
  // ucp_tag_t tag = 0x12;

  // for (int i = 0; i < NOF_RAILS; ++i)
  // {
  //   if (mr_ctx->role == SENDER)
  //   {
  //     ucx_mr_test_send(mr_ctx, i, tag);
  //   }
  //   else
  //   {
  //     ucx_mr_test_recv(mr_ctx, i, tag);
  //   }
  // }
}

ucs_status_t
ucx_mr_create_endpoint(ucx_mr_context_t *mr_ctx)
{
  ucs_status_t status;

  unsigned group_size = mr_ctx->sock_rte_group.size;
  unsigned group_index = mr_ctx->sock_rte_group.is_server ? 0 : 1;
  unsigned peer_index = rte_peer_index(group_size, group_index);

  // printf("group_size %u, group_index %u, peer_index %u\n", group_size, group_index, peer_index);

  int idx = mr_ctx->ep_count;

  if (idx < NOF_RAILS)
  {

    status = ucx_mr_send_local_data(mr_ctx, idx);
    if (status != UCS_OK)
    {
      return status;
    }

    /* receive remote peer's endpoints' data and connect to them */
    status = ucx_mr_recv_remote_data(mr_ctx, peer_index, idx);
    if (status != UCS_OK)
    {
      goto err;
    }

    mr_ctx->ep_count = idx + 1;
  }

  /* sync status across all processes */
  status = ucx_mr_exchange_status(mr_ctx, UCS_OK);
  if (status != UCS_OK)
  {
    goto err_free_eps;
  }

  /* force wireup completion */
  return ucx_mr_flush_workers(mr_ctx);

err_free_eps:
  ucx_mr_destroy_endpoint(mr_ctx);
err:
  (void)ucx_mr_exchange_status(mr_ctx, status);
  return status;
}

ucs_status_t
ucx_mr_create_worker(ucx_mr_context_t *mr_ctx)
{
  ucs_status_t status;
  ucp_worker_params_t *worker_params = (ucp_worker_params_t *)malloc(sizeof(worker_params));
  if (worker_params == NULL)
  {
    ucs_error("failed to allocate memory for worker params");
    status = UCS_ERR_NO_MEMORY;
    return status;
  }

  worker_params->field_mask = UCP_WORKER_PARAM_FIELD_THREAD_MODE;
  worker_params->thread_mode = UCS_THREAD_MODE_SINGLE;

  int idx = mr_ctx->worker_count;

  if (idx < NOF_RAILS)
  {
    status = ucp_worker_create(mr_ctx->ctx[idx], worker_params, &mr_ctx->worker[idx]);
    if (status == UCS_OK)
    {
      mr_ctx->worker_count = idx + 1;
    }
  }
  else
  {
    printf("Maximum number of workers reached!\n");
  }

  free(worker_params);
  return status;
}

ucs_status_t
ucx_mr_init(ucx_mr_context_t *mr_ctx)
{
  ucp_params_t ucp_params;
  ucs_status_t status;
  ucp_config_t *config;

  ucp_params.field_mask = UCP_PARAM_FIELD_FEATURES |
                          UCP_PARAM_FIELD_REQUEST_SIZE |
                          UCP_PARAM_FIELD_REQUEST_INIT |
                          UCP_PARAM_FIELD_NAME;
  ucp_params.features = 0;
  ucp_params.features |= UCP_FEATURE_TAG;
  ucp_params.request_size = sizeof(struct ucx_context);
  ucp_params.request_init = request_init;

int idx = mr_ctx->ctx_count;

  // Set up enviroment and create UCP_Context
  status = ucp_config_read(NULL, NULL, &config);
  if (status != UCS_OK)
  {
    goto err_a;
  }

if (idx==0) {
  status = ucp_config_modify(config, "NET_DEVICES", mr_ctx->rail0);
}
if (idx==1) {
  status = ucp_config_modify(config, "NET_DEVICES", mr_ctx->rail1);
}
if (idx==2) {
  status = ucp_config_modify(config, "NET_DEVICES", mr_ctx->rail2);
}
if (idx==3) {
  status = ucp_config_modify(config, "NET_DEVICES", mr_ctx->rail3);
}

  // if (idx % 2 == 0)
  // {
  //   status = ucp_config_modify(config, "NET_DEVICES", mr_ctx->rail0);
  // }
  // else
  // {
  //   status = ucp_config_modify(config, "NET_DEVICES", mr_ctx->rail1);
  // }
  if (status != UCS_OK)
  {
    goto err_free_config;
  }

  if (idx < NOF_RAILS)
  {
    status = ucp_init(&ucp_params, config, &mr_ctx->ctx[idx]);
    if (status == UCS_OK)
    {
      mr_ctx->ctx_count = idx + 1;
    }
  }
  else
  {
    printf("Maximum number of contexts reached!\n");
  }

err_free_config:
  ucp_config_release(config);
err_a:
  return status;
}

ucs_status_t
ucx_mr_setup(ucx_mr_context_t *mr_ctx)
{
  ucs_status_t status;
  mr_ctx->ctx_count = 0;
  mr_ctx->worker_count = 0;
  mr_ctx->ep_count = 0;

  // Connect Communication
  ucx_mr_connect_comm(mr_ctx);
  DEBUG_PRINT("Connection established!\n");

  // Setup UCX Contexts
  for (int i = 0; i < NOF_RAILS; ++i)
  {
    status = ucx_mr_init(mr_ctx);
    if (status != UCS_OK)
    {
      ucx_mr_cleanup(mr_ctx, CTX);
      return status;
    }
  }
  DEBUG_PRINT("Context created!\n");

  // Setup UCX Worker
  for (int i = 0; i < NOF_RAILS; ++i)
  {
    status = ucx_mr_create_worker(mr_ctx);
    if (status != UCS_OK)
    {
      ucx_mr_cleanup(mr_ctx, WORKER);
      return status;
    }
  }
  DEBUG_PRINT("Worker created!\n");

  // Setup UCX Endpoints
  for (int i = 0; i < NOF_RAILS; ++i)
  {
    status = ucx_mr_create_endpoint(mr_ctx);
    if (status != UCS_OK)
    {
      ucx_mr_cleanup(mr_ctx, EP);
      return status;
    }
  }
  DEBUG_PRINT("Endpoints created!\n");

  // Setup Cuda Peer Access
  cudaSetDevice(0);
  cudaDeviceEnablePeerAccess(1, 0);
  cudaDeviceEnablePeerAccess(2, 0);
  cudaDeviceEnablePeerAccess(3, 0);
  cudaSetDevice(1);
  cudaDeviceEnablePeerAccess(0, 0);
  cudaDeviceEnablePeerAccess(2, 0);
  cudaDeviceEnablePeerAccess(3, 0);
  cudaSetDevice(2);
  cudaDeviceEnablePeerAccess(0, 0);
  cudaDeviceEnablePeerAccess(1, 0);
  cudaDeviceEnablePeerAccess(3, 0);
  cudaSetDevice(3);
  cudaDeviceEnablePeerAccess(0, 0);
  cudaDeviceEnablePeerAccess(1, 0);
  cudaDeviceEnablePeerAccess(2, 0);


  DEBUG_PRINT("Cuda Peer Access enabled!\n");

  int canAccess;

  cudaDeviceCanAccessPeer(&canAccess, 0, 1);
  printf("Can Access 0 -> 1: %d\n", canAccess);

  cudaDeviceCanAccessPeer(&canAccess, 1, 0);
  printf("Can Access 1 -> 0: %d\n", canAccess);

  return status;
}

#endif // UCX_MR_SETUP_H