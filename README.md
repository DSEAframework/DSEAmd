# Data Streaming for Explicit Algorithms - Molecular Dynamics (DSEAmd)

## Introduction
DSEAmd is the implementation of Molecular Dynamics in DSEA.

## Getting Started

### Prerequisites
The following tools are required to build dsea
* cmake 3.21 or newer
* C++ compiler with C++17 support
* CUDA version 12.0 or newer
* MPI library with MPI 2.0 support

Optional requirements are
* UCX sources and libraries version 1.17 or newer for multi rail data transfer between nodes

### Installation

#### makefile
The makefile was used to build the code on the system used for benchmarking. Paths need to be adjusted for other systems.

#### cmak - experimental
To build dsea into a subdirectory `build` execute the following commands from this top level directory

```shell
cmake -B build/ -S .
cmake --build build/
```


### Execution
DSEA works in a ring of processes. When executing the program, care must be taken of the pinning of processes and te assignment of GPUs to processes. Within a node, processes communicate using MPI. Across nodes, the muli rail communication based on UCX can be used. The first and the last process in a node should use multi rail communication via UCX to increase scaling of performance. 

## Citation
Please cite this work as:

"M. Rose, S. Homes, L. Ramsperger, J. Gracia, C. Niethammer, and J. Vrabec. Cyclic Data Streaming on GPUs for Short Range
Stencils Applied to Molecular Dynamics. Submitted to Euro-Par 2025."
