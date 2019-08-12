cmake_minimum_required(VERSION 3.7)

project(clde)
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

set (SOURCES CLmanagerBase.cpp
             CLmanager.cpp
             CLBLAS.cpp)

set (KERNELS kernels/Daxpy.cl)

foreach(kernel ${KERNELS})
    string(REPLACE ".cl" ".h" kernel_header ${kernel})
    list(APPEND KERNEL_HEADERS ${kernel_header})
    add_custom_command(
        OUTPUT  ${kernel_header}
        COMMAND ${PARENT_DIR}/util/clToHeader.py -i ${CMAKE_CURRENT_SOURCE_DIR}/${kernel} -o ${CMAKE_CURRENT_BINARY_DIR}/${kernel_header}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${kernel}
    )
endforeach()

find_package(OpenCL REQUIRED)

include_directories (../include)
include_directories (~/include)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library(${PROJECT_NAME} STATIC ${SOURCES}
                                   ${KERNEL_HEADERS})
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
target_link_libraries(${PROJECT_NAME} OpenCL::OpenCL)
