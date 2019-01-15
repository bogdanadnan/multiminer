cmake_minimum_required (VERSION 2.8)
project (mm_gpu_gate CXX)
set (CMAKE_CXX_STANDARD 11)

set(GPU_GATE_SOURCE
		argon2d-gpu-gate.cpp)

find_package(CUDA)
find_package(OpenCL)

if(CUDA_FOUND AND NOT NO_CUDA)
    include_directories(${CUDA_INCLUDE_DIRS})
	cuda_add_library(mm_gpu_gate SHARED ${GPU_GATE_SOURCE})
    add_definitions(-DHAVE_CUDA=1)
else()
	add_library(mm_gpu_gate SHARED ${GPU_GATE_SOURCE})
endif()

if(OpenCL_FOUND)
    include_directories(${OpenCL_INCLUDE_DIRS})
    add_definitions(-DHAVE_OPENCL=1)
endif()

add_subdirectory(argon2-gpu)

include_directories(${CMAKE_SOURCE_DIR})

target_link_libraries(mm_gpu_gate argon2-cuda argon2-opencl argon2 ${OpenCL_LIBRARIES})
target_compile_definitions(mm_gpu_gate PRIVATE EXPORT_SYMBOLS=1)
target_compile_options(mm_gpu_gate PRIVATE -O3 -mtune=native -march=native -I.)

add_custom_command(
        TARGET mm_gpu_gate POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                "$<TARGET_FILE:mm_gpu_gate>"
                "${CMAKE_BINARY_DIR}/")
