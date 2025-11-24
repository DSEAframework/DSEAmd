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

##### Building on Cray systems with AMD hardware
A separate makefile (makefile_hip_cpe) is provided for building DSEAmd on Cray supercomputers.
The CPE compilers are used. Multi rail communication is currently disabled.
When changing the code, edit the .cu files. All CUDA sources are automatically hipified.

#### cmake - experimental
To build dsea into a sub directory `build` execute the following commands from this top level directory

```shell
cmake -B build/ -S .
cmake --build build/
```

### Execution
DSEA works in a ring of processes. When executing the program, care must be taken of the pinning of processes and the assignment of GPUs to processes.
Within a node, processes communicate using MPI. Across nodes, the muli rail communication based on UCX can be used.
The first and the last process in a node should use multi rail communication via UCX to increase scaling of performance. 

### Visualization
Output is currently generated in vtr format by a call to the routine **caller_output_vtk_rectilinear**.
The generation of output must be enabled in the code.
The output currently provides density of molecules in cells.

## Citation
Please cite this work as:

"M. Rose, S. Homes, L. Ramsperger, J. Gracia, C. Niethammer, and J. Vrabec. Cyclic Data Streaming on GPUs for Short Range
Stencils Applied to Molecular Dynamics. HeteroPar 2025, accepted."

The paper is also available here:
https://arxiv.org/abs/2507.11289