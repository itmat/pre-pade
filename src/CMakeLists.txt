set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_definitions("-DLOG_LEVEL_INFO")

find_library(Z z)
find_library(PT pthread)


# SAMTOOLS
include(FindSamtools)
include_directories(${Samtools_INCLUDE_DIRS})

# add_library(SAM STATIC ${SAM_SOURCES})

set(quantify_SOURCES main.c quant.c)
add_executable(quantify ${quantify_SOURCES})
target_link_libraries(quantify Z)
target_link_libraries(quantify PT)
target_link_libraries(quantify Samtools_LIBRARIES)

set(dumptranscripts_SOURCES quant.c dumptranscripts.c)
add_executable(dumptranscripts ${dumptranscripts_SOURCES})
