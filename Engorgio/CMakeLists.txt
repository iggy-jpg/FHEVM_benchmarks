cmake_minimum_required(VERSION 3.16)
project(ENGORGIO VERSION 1.0 LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 20)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_definitions(-DENGORGIO_DEBUG)
    message(STATUS "START ENGORGIO DEBUG MODE")
    set(CMAKE_CXX_FLAGS "-O0 -g")
else()
    set(CMAKE_CXX_FLAGS "-march=native -O3 -funroll-loops")
endif()


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(WITH_OPENMP OFF CACHE BOOL "Control OpenMP in OpenFHE")


set(OPENFHE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/openfhe-development)


add_subdirectory(${OPENFHE_SOURCE_DIR} ${OPENFHE_SOURCE_DIR}/build)


set(ENGORGIO_SOURCE_FILES "")
add_subdirectory(src/ENGORGIO)

find_package(OpenMP)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif() 

add_library(ENGORGIO STATIC ${ENGORGIO_SOURCE_FILES})
 
target_include_directories(ENGORGIO 
    PUBLIC ${OPENFHE_SOURCE_DIR}/src/core/include
    PUBLIC ${OPENFHE_SOURCE_DIR}/src/pke/include
    PUBLIC ${OPENFHE_SOURCE_DIR}/src/binfhe/include
    PUBLIC ${OPENFHE_SOURCE_DIR}/third-party/cereal/include
    PUBLIC ${OPENFHE_SOURCE_DIR}/build/src/core  
    PUBLIC src
)
 
target_link_directories(ENGORGIO 
    PUBLIC ${CMAKE_BINARY_DIR}/thirdparty/openfhe-development/lib
)
 
target_link_libraries(ENGORGIO 
    PUBLIC OPENFHEcore OPENFHEpke OPENFHEbinfhe
)

message(STATUS "WITH_OPENMP: ${WITH_OPENMP}")
 
add_subdirectory(test)
