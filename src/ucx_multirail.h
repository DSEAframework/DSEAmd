#ifndef UCX_MULTIRAIL_H
#define UCX_MULTIRAIL_H

#include <cuda_runtime.h>
#include <cuda.h>
#include <driver_types.h>

#include <ucp/api/ucp.h>

#include "ucx_mr_types.h"

#include "ucx_mr_aux.h"
#include "ucx_mr_cleanup.h"
#include "ucx_mr_recv.h"
#include "ucx_mr_send.h"
#include "ucx_mr_setup.h"
#include "ucx_mr_sock_comm.h"

#include "quad_simple.h"

#endif // UCX_MULTIRAIL_H