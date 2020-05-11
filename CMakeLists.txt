cmake_minimum_required(VERSION 3.13)

project(powerstrip
        VERSION 0.1
        LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)

# add zstd
set(ZSTD_BUILD_STATIC TRUE)
set(ZSTD_LEGACY_SUPPORT FALSE)
add_subdirectory("extern/zstd/build/cmake")

# add simdcomp
include_directories("${PROJECT_SOURCE_DIR}/extern/simdcomp/include/")
file(GLOB SIMDCOMP_SOURCES "${PROJECT_SOURCE_DIR}/extern/simdcomp/src/*.c")
add_library(simdcomp STATIC ${SIMDCOMP_SOURCES})

include_directories(include)

add_executable(pstrip src/powerstrip.cpp src/compressbin.cpp)
target_link_libraries(pstrip libzstd_static)
target_link_libraries(pstrip simdcomp)

add_executable(unpstrip src/powerstrip.cpp src/decompressbin.cpp)
target_link_libraries(unpstrip libzstd_static)
target_link_libraries(unpstrip simdcomp)

add_library(powerstrip SHARED src/powerstrip.cpp include/powerstrip.h)
target_link_libraries(powerstrip libzstd_static)
target_link_libraries(powerstrip simdcomp)

add_executable(test test/test.cpp)

# add rapidcheck
add_subdirectory("extern/rapidcheck")
target_link_libraries(test rapidcheck)
target_link_libraries(test powerstrip)
