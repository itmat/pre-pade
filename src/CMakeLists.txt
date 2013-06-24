set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_definitions("-DLOG_LEVEL_INFO")

find_library(Z z)
find_library(PT pthread)


# SAMTOOLS
include(FindSamtools)

# SAM/BAM
# set(SAM_DIR "../vendor/samtools")
# file(GLOB SAM_SOURCES ${SAM_DIR}/*.c)
include_directories(${Samtools_INCLUDE_DIRS})

# add_library(SAM STATIC ${SAM_SOURCES})

set(quantify_SOURCES main.c quant.c)
add_executable(quantify ${quantify_SOURCES})

set(dumptranscripts_SOURCES quant.c dumptranscripts.c)
add_executable(dumptranscripts ${dumptranscripts_SOURCES})
