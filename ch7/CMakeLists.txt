cmake_minimum_required(VERSION 2.8)
project(vo1)

#set(CMAKE_BUILD_TYPE "Release")
#add_definitions("-DENABLE_SSE")
#set(CMAKE_CXX_FLAGS "-std=c++11 -O2 ${SSE_FLAGS} -msse4")
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(OpenCV 4 REQUIRED)
Find_Package(G2O REQUIRED)
Find_Package(Eigen3 REQUIRED)
Find_Package(Sophus REQUIRED)
Find_Package(CSparse REQUIRED)

#G2O needs install fmt library:       sudo apt-get install -y libfmt-dev
find_package(fmt)
include_directories(fmt)

SET(G2O_LIBS g2o_csparse_extension g2o_stuff g2o_core cxsparse fmt::fmt)

include_directories(
        ${OpenCV_INCLUDE_DIRS}
        ${G2O_INCLUDE_DIRS}
        ${Sophus_INCLUDE_DIRS}
        "/usr/include/eigen3/"
)

add_executable(orb_cv orb_cv.cpp)
target_link_libraries(orb_cv ${OpenCV_LIBS})

#add_executable(orb_self orb_self.cpp)
#target_link_libraries(orb_self ${OpenCV_LIBS})

# add_executable( pose_estimation_2d2d pose_estimation_2d2d.cpp extra.cpp ) # use this if in OpenCV2 
add_executable(pose_estimation_2d2d pose_estimation_2d2d.cpp)
target_link_libraries(pose_estimation_2d2d ${OpenCV_LIBS})

# # add_executable( triangulation triangulation.cpp extra.cpp) # use this if in opencv2
add_executable(triangulation triangulation.cpp)
target_link_libraries(triangulation ${OpenCV_LIBS})

add_executable(pose_estimation_3d2d pose_estimation_3d2d.cpp)
target_link_libraries(pose_estimation_3d2d
        ${G2O_LIBS}
        ${OpenCV_LIBS})

add_executable(pose_estimation_3d3d pose_estimation_3d3d.cpp)
target_link_libraries(pose_estimation_3d3d
        ${G2O_LIBS}
        ${OpenCV_LIBS})
