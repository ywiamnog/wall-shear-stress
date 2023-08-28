project(deformation)
cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 17)

#set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_BUILD_TYPE Release)
#set (CMAKE_CXX_FLAGS "-Wall")
#set (CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wunused")

add_executable(deformation "source/main.cpp")

target_include_directories(deformation PUBLIC "include/")
target_include_directories(deformation PUBLIC "../../core_tools/voxels/include/")
#target_include_directories(undo_axis_alignment PUBLIC "/home/yimo/git-repos/json/single_include/")

add_subdirectory(../../core_tools build_core_tools)
target_include_directories(deformation PUBLIC "../../core_tools/common/include")
target_include_directories(deformation PUBLIC "../../core_tools/affinities/include")
target_include_directories(deformation PUBLIC "../../core_tools/misc_io/include")

target_link_libraries(deformation PUBLIC common)
target_link_libraries(deformation PUBLIC misc_io)


