cmake_minimum_required(VERSION 3.5)

project(argon2-gpu CXX)
set(BINARY_INSTALL_DIR /usr/local/bin)
set(LIBRARY_INSTALL_DIR /usr/local/lib)
set(INCLUDE_INSTALL_DIR /usr/local/include)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()

set(CUDA_FOUND FALSE)

if(NOT NO_CUDA)
    find_package(CUDA)
endif()

if(CUDA_FOUND)
    message("INFO: Using CUDA version ${CUDA_VERSION}")
    add_definitions(-DHAVE_CUDA=1)
else()
    message("INFO: Building without CUDA support")
    add_definitions(-DHAVE_CUDA=0)
endif()

if(CUDA_FOUND)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};
        -std=c++11;-O3;--ptxas-options=-v;-arch sm_30;-lineinfo
    )
endif()

add_subdirectory(ext/argon2)

add_library(argon2-gpu-common SHARED
    lib/argon2-gpu-common/argon2params.cpp
    lib/argon2-gpu-common/blake2b.cpp
)
target_include_directories(argon2-gpu-common INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_include_directories(argon2-gpu-common PRIVATE
    include/argon2-gpu-common
    lib/argon2-gpu-common
)

include_directories(include)

if(CUDA_FOUND)
    cuda_add_library(argon2-cuda SHARED
        lib/argon2-cuda/device.cpp
        lib/argon2-cuda/globalcontext.cpp
        lib/argon2-cuda/programcontext.cpp
        lib/argon2-cuda/processingunit.cpp
            lib/argon2-cuda/kernels.cu
    )
    include_directories(${CUDA_INCLUDE_DIRS})
else()
    add_library(argon2-cuda SHARED
        lib/argon2-cuda/nocuda.cpp)
endif()

target_include_directories(argon2-cuda PRIVATE
    include/argon2-cuda
    lib/argon2-cuda
)
target_include_directories(argon2-cuda INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(argon2-cuda argon2-gpu-common)

add_library(argon2-opencl SHARED
    lib/argon2-opencl/device.cpp
    lib/argon2-opencl/globalcontext.cpp
    lib/argon2-opencl/kernelloader.cpp
    lib/argon2-opencl/programcontext.cpp
    lib/argon2-opencl/processingunit.cpp
    lib/argon2-opencl/kernelrunner.cpp
)
target_include_directories(argon2-opencl INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_include_directories(argon2-opencl PRIVATE
    include/argon2-opencl
    lib/argon2-opencl
)
target_link_libraries(argon2-opencl
    argon2-gpu-common -lOpenCL
)

install(
    TARGETS argon2-gpu-common argon2-opencl argon2-cuda
    DESTINATION ${LIBRARY_INSTALL_DIR}
)
install(FILES
    include/argon2-gpu-common/argon2-common.h
    include/argon2-gpu-common/argon2params.h
    include/argon2-opencl/cl.hpp
    include/argon2-opencl/opencl.h
    include/argon2-opencl/device.h
    include/argon2-opencl/globalcontext.h
    include/argon2-opencl/programcontext.h
    include/argon2-opencl/processingunit.h
    include/argon2-opencl/kernelrunner.h
    include/argon2-cuda/cudaexception.h
    include/argon2-cuda/kernels.h
    include/argon2-cuda/device.h
    include/argon2-cuda/globalcontext.h
    include/argon2-cuda/programcontext.h
    include/argon2-cuda/processingunit.h
    DESTINATION ${INCLUDE_INSTALL_DIR}
)
