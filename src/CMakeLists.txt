cmake_minimum_required(VERSION 3.7)

project(clde)
set (SOURCES CLmanagerBase.cpp)
find_package(OpenCL REQUIRED)

include_directories (../include)


add_library(${PROJECT_NAME} STATIC ${SOURCES})
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
target_link_libraries(${PROJECT_NAME} OpenCL::OpenCL)
