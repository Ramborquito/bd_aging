cmake_minimum_required(VERSION 3.31)
project(bd_aging CUDA)

set(CMAKE_CUDA_STANDARD 20)

add_executable(bd_aging.bin main.cu
        cell_locate.cu
        colors.cu
        Contracting_Vector3D.cu
        get_forces.cu
        get_gder.cu
        potentials.cu
        rannew64.cu
        update_ermak.cu
        set_temp.cu)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

set_target_properties(bd_aging.bin PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(bd_aging.bin ${GSL_LIBRARIES})
