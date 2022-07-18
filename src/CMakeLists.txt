set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(BINARY_INCLUDE "${PROJECT_BINARY_DIR}/include")
set(SOURCE_FILES main.cpp)
set(INCLUDE_DIRECTORIES "${BINARY_INCLUDE}")


# Add Eigen3:
set(HIGHFIVE_USE_EIGEN ON)
#message("CMAKE_SYSTEM_PREFIX_PATH=${CMAKE_SYSTEM_PREFIX_PATH}")
find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "Found my Eigen ${Eigen3_VERSION}: ${EIGEN3_INCLUDE_DIRS}")
else()
    find_package(PkgConfig)
    pkg_check_modules(EIGEN3 REQUIRED eigen3)
endif()

include(FetchContent)
# Add High Five:
find_package(HDF5)
FetchContent_Declare(
        HighFive
        GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
        GIT_TAG 3c8796ff1c849a6056ebd65fe540e4299799f09e
)
FetchContent_MakeAvailable(HighFive)
set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} "${highfive_SOURCE_DIR}/include")

# Set the GPU to compile for
set(GPU_TARGETS "gfx1012" CACHE STRING "GPU targets to compile for")

# Search for rocm in common locations
list(APPEND CMAKE_PREFIX_PATH /opt/rocm/hip /opt/rocm)
# Find hip
# find_package(hip)
if(hip_FOUND)
    message(STATUS "Found my hip ${hip_VERSION}: ${hip_INCLUDE_DIRS}")
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} "${hip_INCLUDE_DIRS}")
endif()

#find openmp
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found my OpenMP ${OpenMP_CXX_VERSION}: ${OpenMP_CXX_INCLUDE_DIRS}  ${OpenMP_CXX_LIBRARIES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")




# configure a header file to pass some of the CMake settings
# to the source code
configure_file("SrcConfig.h.in" "${BINARY_INCLUDE}/SrcConfig.h")

add_subdirectory(LeapFrogIntegrator)
add_subdirectory(IsingModel)
add_subdirectory(HMCGenerator)
add_subdirectory(MyTypes)
add_subdirectory(MultiLevelHMCGenerator)
add_subdirectory(BaseModel)
add_subdirectory(Analyzer)
add_subdirectory(XYModel)


# add the executable
add_executable("${PROJECT_NAME}" ${SOURCE_FILES})

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
target_include_directories("${PROJECT_NAME}" PUBLIC ${INCLUDE_DIRECTORIES})
target_link_libraries("${PROJECT_NAME}" PUBLIC ${HDF5_LIBRARIES} Eigen3::Eigen OpenMP::OpenMP_CXX)
