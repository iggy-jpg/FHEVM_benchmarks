add_executable(comparison_test comparison_test.cpp)
target_include_directories(comparison_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(comparison_test ENGORGIO)

add_executable(sort_test sort_test.cpp)
target_include_directories(sort_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(sort_test ENGORGIO)

add_executable(hybrid_query_test_1 hybrid_query_test_1.cpp)
target_include_directories(hybrid_query_test_1 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(hybrid_query_test_1 ENGORGIO)

add_executable(hybrid_query_test_2 hybrid_query_test_2.cpp)
target_include_directories(hybrid_query_test_2 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(hybrid_query_test_2 ENGORGIO)

add_executable(relational_query_test relational_query_test.cpp)
target_include_directories(relational_query_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(relational_query_test ENGORGIO)

add_executable(vectorized_query_test_1 vectorized_query_test_1.cpp)
target_include_directories(vectorized_query_test_1 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(vectorized_query_test_1 ENGORGIO)

add_executable(vectorized_query_test_2 vectorized_query_test_2.cpp)
target_include_directories(vectorized_query_test_2 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(vectorized_query_test_2 ENGORGIO)

add_executable(topk_test topk_test.cpp)
target_include_directories(topk_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(topk_test ENGORGIO)

add_executable(sync_test sync_test.cpp)
target_include_directories(sync_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(sync_test ENGORGIO)

add_executable(minimax_test compsigneval.cpp)
target_include_directories(minimax_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(minimax_test ENGORGIO)

add_executable(transaction_test hybrid_transaction_test.cpp)
target_include_directories(transaction_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(transaction_test ENGORGIO)

add_executable(operation_performance_test operation_performance_test.cpp)
target_include_directories(operation_performance_test
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(operation_performance_test ENGORGIO)
