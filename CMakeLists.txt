cmake_minimum_required(VERSION 3.22)

# set the project name and version
project(MultiLevelHMC VERSION 1.0)

set(SOURCE_DIR_NAME "src")
set(DOCUMENTATION_DIR_NAME "docs")
set(DATA_DIR_NAME "data")

# Print some info:
message("PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
message("PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}")


# Options:
option(BUILD_DOC "Build documentation" ON)

# Compile src
add_subdirectory(${SOURCE_DIR_NAME})

# Docs
find_package(Doxygen)
if (BUILD_DOC)
    if (DOXYGEN_FOUND)
        add_subdirectory(${DOCUMENTATION_DIR_NAME})
    else ()
        message("Doxygen needs to be installed to build the documentation!")
    endif ()
endif ()




