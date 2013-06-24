cmake_minimum_required (VERSION 2.6)
set(CMAKE_C_COMPILER gcc)

# Put the executables in the bin dir.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project (PREPADE C)


# force an out-of-source build of the code
set(CMAKE_MODULE_PATH ${PREPADE_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})
include(MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("Please build ${PROJECT_NAME} out-of-source")

# includes
include_directories("${PROJECT_SOURCE_DIR}/include")
add_subdirectory(src)

