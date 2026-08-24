cmake_minimum_required(VERSION 3.16)
project(HomGateProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# AppleClang on arm64 rejects -march=native; -mcpu=native is the equivalent there.
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    set(AUTOHOG_ARCH_FLAG "-mcpu=native")
else()
    set(AUTOHOG_ARCH_FLAG "-march=native")
endif()
set(CMAKE_CXX_FLAGS "${AUTOHOG_ARCH_FLAG} -O3 -g -funroll-loops -Wall -Wextra -pedantic -Wno-sign-compare")

# OpenMP (on macOS/AppleClang, hint with -DOpenMP_ROOT=$(brew --prefix libomp))
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    if(OpenMP_CXX_INCLUDE_DIRS)
        include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    endif()
    link_libraries(${OpenMP_CXX_LIBRARIES})
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(AUTOHOG_ZK_BACKEND "blst" CACHE STRING "ZK elliptic-curve backend: mcl or blst")
set_property(CACHE AUTOHOG_ZK_BACKEND PROPERTY STRINGS mcl blst)
if(NOT AUTOHOG_ZK_BACKEND STREQUAL "mcl" AND NOT AUTOHOG_ZK_BACKEND STREQUAL "blst")
    message(FATAL_ERROR "AUTOHOG_ZK_BACKEND must be either 'mcl' or 'blst'")
endif()

set(AUTOHOG_ZK_CURVE "bls12_381" CACHE STRING "ZK elliptic-curve: bls12_381 or bn_snark1")
set_property(CACHE AUTOHOG_ZK_CURVE PROPERTY STRINGS bls12_381 bn_snark1)
if(NOT AUTOHOG_ZK_CURVE STREQUAL "bls12_381" AND NOT AUTOHOG_ZK_CURVE STREQUAL "bn_snark1")
    message(FATAL_ERROR "AUTOHOG_ZK_CURVE must be either 'bls12_381' or 'bn_snark1'")
endif()

if(AUTOHOG_ZK_BACKEND STREQUAL "blst" AND AUTOHOG_ZK_CURVE STREQUAL "bn_snark1")
    message(FATAL_ERROR "BLST backend does not support bn_snark1 curve. Use bls12_381 or switch to mcl backend.")
endif()

message(STATUS "ZK Backend: ${AUTOHOG_ZK_BACKEND}")
message(STATUS "ZK Curve: ${AUTOHOG_ZK_CURVE}")

# Homebrew default prefix is not on the compiler search path on macOS arm64
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    include_directories(SYSTEM /opt/homebrew/include)
    link_directories(/opt/homebrew/lib)
endif()

# Add TFHEpp library
add_subdirectory(thirdparties/TFHEpp)

# TFHEpp exports its compile definitions (FFT backend selection etc.) via
# TFHEpp_DEFINITIONS; consumers must apply them or headers pick the wrong backend.
if(TFHEpp_DEFINITIONS)
    add_compile_definitions(${TFHEpp_DEFINITIONS})
endif()

# Include directories
include_directories(thirdparties/TFHEpp/src)
include_directories(include)

if(AUTOHOG_ZK_BACKEND STREQUAL "mcl")
    add_subdirectory(thirdparties/mcl)
    set(AUTOHOG_ZK_BACKEND_LIB mcl::mcl)
elseif(AUTOHOG_ZK_BACKEND STREQUAL "blst")
    set(BLST_ROOT "${PROJECT_SOURCE_DIR}/thirdparties/blst" CACHE PATH
        "Path to a blst source tree")
    if(NOT EXISTS "${BLST_ROOT}/bindings")
        message(FATAL_ERROR "BLST_ROOT does not contain bindings/: ${BLST_ROOT}")
    endif()
    set(BLST_LIBRARY "${BLST_ROOT}/libblst.a" CACHE FILEPATH
        "Path to a prebuilt blst static library")
    option(AUTOHOG_USE_PREBUILT_BLST
        "Use BLST_LIBRARY as-is instead of invoking blst/build.sh" OFF)

    if(AUTOHOG_USE_PREBUILT_BLST)
        if(NOT EXISTS "${BLST_LIBRARY}")
            message(FATAL_ERROR
                "AUTOHOG_USE_PREBUILT_BLST=ON but library is missing: ${BLST_LIBRARY}")
        endif()
    else()
        add_custom_command(
            OUTPUT ${BLST_LIBRARY}
            COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} sh build.sh
            WORKING_DIRECTORY ${BLST_ROOT}
            DEPENDS
                ${BLST_ROOT}/build.sh
                ${BLST_ROOT}/src/server.c
                ${BLST_ROOT}/build/assembly.S
            COMMENT "Building blst"
            VERBATIM
        )
        add_custom_target(blst-build DEPENDS ${BLST_LIBRARY})
    endif()
    add_library(blst STATIC IMPORTED GLOBAL)
    set_target_properties(blst PROPERTIES
        IMPORTED_LOCATION ${BLST_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES ${BLST_ROOT}/bindings
    )
    if(NOT AUTOHOG_USE_PREBUILT_BLST)
        add_dependencies(blst blst-build)
    endif()
    set(AUTOHOG_ZK_BACKEND_LIB blst)
endif()




find_package(OpenSSL REQUIRED)


file(GLOB HOMOPS_SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")


add_library(homops STATIC ${HOMOPS_SRCS})

target_include_directories(homops PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/thirdparties/TFHEpp/include
    ${PROJECT_SOURCE_DIR}/thirdparties/TFHEpp/src
    ${PROJECT_SOURCE_DIR}/thirdparties/TFHEpp/spqlios/spqlios
)

target_compile_definitions(homops PUBLIC AUTOHOG_ZK_BACKEND_${AUTOHOG_ZK_BACKEND}=1)
target_compile_definitions(homops PUBLIC AUTOHOG_ZK_CURVE_${AUTOHOG_ZK_CURVE}=1)
target_link_libraries(homops PUBLIC tfhe++ ${AUTOHOG_ZK_BACKEND_LIB}
                                   OpenSSL::SSL OpenSSL::Crypto)

add_library(fhe-api SHARED src/fhe_key_api.cpp)
target_link_libraries(fhe-api PUBLIC homops)

# Find and link OpenMP
find_package(OpenMP)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(Coprocessor-operator)
