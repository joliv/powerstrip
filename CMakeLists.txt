cmake_minimum_required(VERSION 3.16)
project(powerstrip
        VERSION 0.1
        LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)

# add zstd
set(ZSTD_BUILD_STATIC TRUE)
set(ZSTD_LEGACY_SUPPORT FALSE)
add_subdirectory("extern/zstd/build/cmake")

# add simdcomp
# TODO get libsimdcomp.a out of the root dir
include_directories("${PROJECT_SOURCE_DIR}/extern/simdcomp/include/")
file(GLOB SIMDCOMP_SOURCES "${PROJECT_SOURCE_DIR}/extern/simdcomp/src/*.c")
add_library(simdcomp STATIC ${SIMDCOMP_SOURCES})

add_executable(compress src/compress.cpp)
target_link_libraries(compress libzstd_static)
target_link_libraries(compress simdcomp)

# compile-time options, if needed
set(DEBUG_OPTIONS, "-Wall")
set(RELEASE_OPTIONS, "-Wall")
target_compile_options(compress PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(compress PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")

add_executable(decompress src/decompress.cpp)
