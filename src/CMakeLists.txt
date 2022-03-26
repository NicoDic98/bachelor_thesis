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
        GIT_TAG 5db7df453cd32811cd63b75627236cf13d337b66
)
FetchContent_MakeAvailable(HighFive)
set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} "${highfive_SOURCE_DIR}/include")


# configure a header file to pass some of the CMake settings
# to the source code
configure_file("SrcConfig.h.in" "${BINARY_INCLUDE}/SrcConfig.h")

add_subdirectory(LeapFrogIntegrator)
add_subdirectory(IsingModel)
add_subdirectory(HMCGenerator)
add_subdirectory(MyTypes)
add_subdirectory(MultiLevelHMCGenerator)


# add the executable
add_executable("${PROJECT_NAME}" ${SOURCE_FILES})

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
target_include_directories("${PROJECT_NAME}" PUBLIC ${INCLUDE_DIRECTORIES})
target_link_libraries("${PROJECT_NAME}" PUBLIC ${HDF5_LIBRARIES} Eigen3::Eigen)
