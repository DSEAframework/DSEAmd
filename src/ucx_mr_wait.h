#ifndef UCX_MR_WAIT_H
#define UCX_MR_WAIT_H

#include <ucp/api/ucp.h>

static ucs_status_t ucx_wait(ucp_worker_h ucp_worker, struct ucx_context *request, const char *op_str, const char *data_str)
{
    ucs_status_t status;
    if (UCS_PTR_IS_ERR(request))
    {
        status = UCS_PTR_STATUS(request);
    }
    else if (UCS_PTR_IS_PTR(request))
    {
        while (!request->completed)
        {
            ucp_worker_progress(ucp_worker);
        }

        request->completed = 0;
        status = ucp_request_check_status(request);
        ucp_request_free(request);
    }
    else
    {
        status = UCS_OK;
    }

    if (status != UCS_OK)
    {
        fprintf(stderr, "unable to %s %s (%s)\n", op_str, data_str,
                ucs_status_string(status));
    }
    else
    {
        DEBUG_PRINT("finish to %s %s\n", op_str, data_str);
    }

    return status;
}


static ucs_status_t ucx_wait_max(ucp_worker_h ucp_worker, struct ucx_context *request, const char *op_str, const char *data_str, int n_progress_max) {
	ucs_status_t status;
	if (UCS_PTR_IS_ERR(request)) {
		status = UCS_PTR_STATUS(request);
	}
	else if (UCS_PTR_IS_PTR(request)) {
		int n_progress=0;
		while ((!request->completed)&&(n_progress<n_progress_max)) {
			ucp_worker_progress(ucp_worker);
			n_progress++;
		}

		if (n_progress==n_progress_max) {
			// max nmber of tries exceede, the operation might have completed in the last call of ucp_worker_progress
			status = ucp_request_check_status(request);
			if (status == UCS_OK) {
				request->completed = 0;
				ucp_request_free(request);
			}
			return status;
		}

		request->completed = 0;
		status = ucp_request_check_status(request);
		ucp_request_free(request);
	}
	else {
		status = UCS_OK;
	}

	if (status != UCS_OK) {
		fprintf(stderr, "unable to %s %s (%s)\n", op_str, data_str,
				ucs_status_string(status));
	}
	else {
		DEBUG_PRINT("finish to %s %s\n", op_str, data_str);
	}

	return status;
}

#endif // UCX_MR_WAIT_H