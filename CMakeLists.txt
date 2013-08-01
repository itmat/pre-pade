cmake_minimum_required (VERSION 2.8)
set(CMAKE_C_COMPILER gcc)
project (PREPADE C)

# Put the executables in the bin dir.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# includes
include_directories("${PROJECT_SOURCE_DIR}/include")
add_subdirectory(src)

# force an out-of-source build of the code
set(CMAKE_MODULE_PATH ${PREPADE_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})
include(MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("Please build ${PROJECT_NAME} out-of-source")

enable_testing()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/testdata/arabidopsis.gtf
    ${CMAKE_BINARY_DIR}/testdata/arabidopsis.gtf COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/testdata/RUM.sam
    ${CMAKE_BINARY_DIR}/testdata/RUM.sam COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/testdata/cigar_bug.sam
    ${CMAKE_BINARY_DIR}/testdata/cigar_bug.sam COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/testdata/test_next_fragment_paired.sam
    ${CMAKE_BINARY_DIR}/testdata/test_next_fragment_paired.sam COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/testdata/test_next_fragment_single.sam
    ${CMAKE_BINARY_DIR}/testdata/test_next_fragment_single.sam COPYONLY)


add_test(GeneIndex ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testgeneindex)
add_test(SamUtils ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testsamutils)
