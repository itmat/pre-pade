set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_definitions("-DLOG_LEVEL_INFO")

find_library(Z z)
find_library(PT pthread)


# SAMTOOLS
# include(FindSamtools)
# include_directories(Samtools_INCLUDE_DIR)
set(SAMTOOLS_DIR ${PROJECT_SOURCE_DIR}/vendor/samtools)
set(SAMTOOLS_LIBRARY ${SAMTOOLS_DIR}/libbam.a)
add_custom_target(build_samtools ALL
                    COMMAND ${CMAKE_MAKE_PROGRAM} lib
                    WORKING_DIRECTORY ${SAMTOOLS_DIR}
                    COMMENT "Make static compiled samtools library"
                )
add_library(samtools STATIC IMPORTED)
set_property(TARGET samtools APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(samtools PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${SAMTOOLS_LIBRARY}")
add_dependencies(samtools build_samtools)

include_directories(${SAMTOOLS_DIR})

set(quantify_SOURCES quant.c main.c)
add_executable(quantify ${quantify_SOURCES})
target_link_libraries(quantify ${Z})
target_link_libraries(quantify ${PT})
target_link_libraries(quantify samtools)
# target_link_libraries(quantify Samtools_LIBRARIES)
# target_include_directories(quantify INTERFACE ${Samtools_INCLUDE_DIR})
set(dumptranscripts_SOURCES quant.c dumptranscripts.c)
add_executable(dumptranscripts ${dumptranscripts_SOURCES})
target_link_libraries(dumptranscripts ${Z})
target_link_libraries(dumptranscripts ${PT})
target_link_libraries(dumptranscripts samtools)
