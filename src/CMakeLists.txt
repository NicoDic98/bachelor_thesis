set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(BINARY_INCLUDE "${PROJECT_BINARY_DIR}/include")
set(SOURCE_FILES main.cpp)
set(INCLUDE_DIRECTORIES "${BINARY_INCLUDE}")


# Add High Five:
find_package(HDF5)
include(FetchContent)
FetchContent_Declare(
        HighFive
        GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
        GIT_TAG 5db7df453cd32811cd63b75627236cf13d337b66
)
FetchContent_MakeAvailable(HighFive)
message("highfive_SOURCE_DIR=${highfive_SOURCE_DIR}")
message("highfive_POPULATED=${highfive_POPULATED}")
message("highfive_BINARY_DIR=${highfive_BINARY_DIR}")
set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} "${highfive_SOURCE_DIR}/include")


# configure a header file to pass some of the CMake settings
# to the source code
configure_file("SrcConfig.h.in" "${BINARY_INCLUDE}/SrcConfig.h")

add_subdirectory(LeapFrogIntegrator)
add_subdirectory(IsingModel)
add_subdirectory(HMCGenerator)
add_subdirectory(MyTypes)

message("SOURCE_FILES=${SOURCE_FILES}")
message("INCLUDE_DIRECTORIES=${INCLUDE_DIRECTORIES}")

# add the executable
add_executable("${PROJECT_NAME}" ${SOURCE_FILES})

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
target_include_directories("${PROJECT_NAME}" PUBLIC ${INCLUDE_DIRECTORIES})
target_link_libraries("${PROJECT_NAME}" PUBLIC ${HDF5_LIBRARIES})
