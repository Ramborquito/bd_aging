cmake_minimum_required(VERSION 3.31)
project(bd_aging CUDA)

set(CMAKE_CUDA_STANDARD 20)

add_executable(bd_aging_self.bin main.cu
        cell_locate.cu
        get_forces_poly.cu
        get_gder_poly.cu
        potentials.cu
        rannew64.cu
        update_ermak.cu
        fself.cu
        set_temp.cu)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

set_target_properties(bd_aging_self.bin PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(bd_aging_self.bin ${GSL_LIBRARIES})
