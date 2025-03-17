cmake_minimum_required(VERSION 3.21)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.30)
    cmake_policy(SET CMP0167 NEW)
endif()

# Minimum CUDA architecture is 61 because of use of atomicAdd(double *, double)
set(CMAKE_CUDA_ARCHITECTURES "61")

project(dsea
        VERSION 0.1.0
        DESCRIPTION "Data Streaming for Explicit Algorithms - DSEA"
        LANGUAGES C CUDA CXX
        )

find_package(Boost COMPONENTS filesystem REQUIRED )
find_package(CUDAToolkit REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(UCX 1.17 REQUIRED)


add_executable(dsea)

target_sources(dsea
    PRIVATE
        main.cpp
        dsea.cpp
        dsea.cpp
        dsea_block.cpp
        dsea_gpu.cpp
        dsea_input.cpp
        dsea_output.cpp
        dsea_storage.cpp
        dsea_main.cpp
        dsea_kernel.cu
    )
target_link_libraries(dsea
    PUBLIC
        MPI::MPI_CXX
        CUDA::cudart
        OpenMP::OpenMP_CXX
        Boost::filesystem
    )


add_executable(dsea_ucx_send)
target_compile_definitions(dsea_ucx_send PRIVATE MRUCX_SEND )
target_sources(dsea_ucx_send
    PRIVATE
        main.cpp
        dsea.cpp
        dsea.cpp
        dsea_block.cpp
        dsea_gpu.cpp
        dsea_input.cpp
        dsea_output_ucx.cpp
        dsea_storage.cpp
        dsea_main.cpp
        dsea_kernel.cu
    )
target_link_libraries(dsea_ucx_send
    PUBLIC
        MPI::MPI_CXX
        CUDA::cudart
        OpenMP::OpenMP_CXX
        Boost::filesystem
        ucx::ucp
        ucx::ucs
    )


add_executable(dsea_ucx_recv)
target_compile_definitions(dsea_ucx_recv PRIVATE MRUCX_REC )
target_sources(dsea_ucx_recv
    PRIVATE
        main.cpp
        dsea.cpp
        dsea.cpp
        dsea_block.cpp
        dsea_gpu.cpp
        dsea_input_ucx.cpp
        dsea_output.cpp
        dsea_storage.cpp
        dsea_main.cpp
        dsea_kernel.cu
    )
target_link_libraries(dsea_ucx_recv
    PUBLIC
        MPI::MPI_CXX
        CUDA::cudart
        OpenMP::OpenMP_CXX
        Boost::filesystem
        ucx::ucp
        ucx::ucs
    )
