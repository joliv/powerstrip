cmake_minimum_required(VERSION 3.13)
#cmake_minimum_required(VERSION 3.13)
project(powerstrip
        VERSION 0.1
        LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_VERBOSE_MAKEFILE on)

# create lib, bin dirs
#add_custom_command(OUTPUT lib
#        PRE_BUILD
#        COMMAND ${CMAKE_COMMAND} -E make_directory lib)
##add_custom_command(OUTPUT bin
##        PRE_BUILD
##        COMMAND ${CMAKE_COMMAND} -E make_directory bin)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
##set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PROJECT_SOURCE_DIR}/bin/debug)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}/bin/release)

# add zstd
set(ZSTD_BUILD_STATIC TRUE)
set(ZSTD_LEGACY_SUPPORT FALSE)
add_subdirectory("extern/zstd/build/cmake")

# add simdcomp
# TODO get libsimdcomp.a out of the root dir
include_directories("${PROJECT_SOURCE_DIR}/extern/simdcomp/include/")
file(GLOB SIMDCOMP_SOURCES "${PROJECT_SOURCE_DIR}/extern/simdcomp/src/*.c")
add_library(simdcomp STATIC ${SIMDCOMP_SOURCES})

include_directories(include)

add_executable(compress src/powerstrip.cpp src/compressbin.cpp)
target_link_libraries(compress libzstd_static)
target_link_libraries(compress simdcomp)

add_executable(decompress src/powerstrip.cpp src/decompressbin.cpp)
target_link_libraries(decompress libzstd_static)
target_link_libraries(decompress simdcomp)

add_library(powerstrip SHARED src/powerstrip.cpp include/powerstrip.h)
target_link_libraries(powerstrip libzstd_static)
target_link_libraries(powerstrip simdcomp)

add_executable(test test/test.cpp)

# add rapidcheck
add_subdirectory("extern/rapidcheck")
target_link_libraries(test rapidcheck)
target_link_libraries(test powerstrip)

# compile-time options, if needed
#set(DEBUG_OPTIONS, "-Wall")
#set(RELEASE_OPTIONS, "-Wall -O3")
#target_compile_options(compress PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
#target_compile_options(compress PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")
#target_compile_options(decompress PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
#target_compile_options(decompress PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")
