cmake_minimum_required(VERSION 3.16)

project(ks256_cpu LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenFHE CONFIG REQUIRED)

message(STATUS "OpenFHE version: ${BASE_OPENFHE_VERSION}")
message(STATUS "OpenFHE include: ${OpenFHE_INCLUDE}")
message(STATUS "OpenFHE libraries: ${OpenFHE_SHARED_LIBRARIES}")

add_executable(tpch_q6_16bit_cpu
    src/tpch_q6_16bit_cpu.cpp
    src/tpch_boolean_greater.cpp
    src/tpch_boolean_and.cpp
    src/tpch_boolean_multiplier.cpp
    src/tpch_boolean_adder.cpp
)

add_executable(boolean_operator_benchmark_cpu
    src/boolean_operator_benchmark_cpu.cpp
)

add_executable(transaction_64bit_cpu
    src/transaction_64bit_cpu.cpp
)

target_include_directories(tpch_q6_16bit_cpu
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(boolean_operator_benchmark_cpu
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(transaction_64bit_cpu
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(tpch_q6_16bit_cpu
    SYSTEM PRIVATE
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
    ${OpenFHE_INCLUDE}/third-party/include
)

target_include_directories(boolean_operator_benchmark_cpu
    SYSTEM PRIVATE
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
    ${OpenFHE_INCLUDE}/third-party/include
)

target_include_directories(transaction_64bit_cpu
    SYSTEM PRIVATE
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
    ${OpenFHE_INCLUDE}/third-party/include
)

target_link_directories(tpch_q6_16bit_cpu
    PRIVATE
    ${OpenFHE_LIBDIR}
)

target_link_directories(boolean_operator_benchmark_cpu
    PRIVATE
    ${OpenFHE_LIBDIR}
)

target_link_directories(transaction_64bit_cpu
    PRIVATE
    ${OpenFHE_LIBDIR}
)

target_link_libraries(tpch_q6_16bit_cpu
    PRIVATE
    ${OpenFHE_SHARED_LIBRARIES}
)

target_link_libraries(boolean_operator_benchmark_cpu
    PRIVATE
    ${OpenFHE_SHARED_LIBRARIES}
)

target_link_libraries(transaction_64bit_cpu
    PRIVATE
    ${OpenFHE_SHARED_LIBRARIES}
)

separate_arguments(OPENFHE_COMPILE_OPTIONS NATIVE_COMMAND "${OpenFHE_CXX_FLAGS}")
separate_arguments(OPENFHE_LINK_OPTIONS NATIVE_COMMAND "${OpenFHE_EXE_LINKER_FLAGS}")
target_compile_options(tpch_q6_16bit_cpu PRIVATE ${OPENFHE_COMPILE_OPTIONS})
target_link_options(tpch_q6_16bit_cpu PRIVATE ${OPENFHE_LINK_OPTIONS})
target_compile_options(boolean_operator_benchmark_cpu PRIVATE ${OPENFHE_COMPILE_OPTIONS})
target_link_options(boolean_operator_benchmark_cpu PRIVATE ${OPENFHE_LINK_OPTIONS})
target_compile_options(transaction_64bit_cpu PRIVATE ${OPENFHE_COMPILE_OPTIONS})
target_link_options(transaction_64bit_cpu PRIVATE ${OPENFHE_LINK_OPTIONS})

if(UNIX AND NOT APPLE)
    target_link_options(tpch_q6_16bit_cpu PRIVATE -Wl,--no-as-needed)
    target_link_options(boolean_operator_benchmark_cpu PRIVATE -Wl,--no-as-needed)
    target_link_options(transaction_64bit_cpu PRIVATE -Wl,--no-as-needed)
endif()

if(MSVC)
    target_compile_options(tpch_q6_16bit_cpu PRIVATE /W4)
    target_compile_options(boolean_operator_benchmark_cpu PRIVATE /W4)
    target_compile_options(transaction_64bit_cpu PRIVATE /W4)
else()
    target_compile_options(tpch_q6_16bit_cpu PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(boolean_operator_benchmark_cpu PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(transaction_64bit_cpu PRIVATE -Wall -Wextra -Wpedantic)
endif()
