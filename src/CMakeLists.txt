set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(BINARY_INCLUDE "${PROJECT_BINARY_DIR}/include")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file("SrcConfig.h.in" "${BINARY_INCLUDE}/SrcConfig.h")

set(SOURCE_FILES main.cpp)

add_subdirectory(LeapFrogIntegrator)
add_subdirectory(IsingModel)

message("SOURCE_FILES=${SOURCE_FILES}")

# add the executable
add_executable("${PROJECT_NAME}" ${SOURCE_FILES})

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
target_include_directories("${PROJECT_NAME}" PUBLIC
        "${BINARY_INCLUDE}"
        )
