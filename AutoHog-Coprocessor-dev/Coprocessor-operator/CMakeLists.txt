file(GLOB CPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(file_path IN LISTS CPP_FILES)
    get_filename_component(file_name ${file_path} NAME_WE)
    add_executable(${file_name} ${file_path})

    target_include_directories(${file_name} PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/thirdparties/TFHEpp/include
        ${PROJECT_SOURCE_DIR}/thirdparties/TFHEpp/src
        ${PROJECT_SOURCE_DIR}/thirdparties/TFHEpp/spqlios/spqlios
    )

    target_link_libraries(${file_name} PRIVATE homops tfhe++)
    target_link_libraries(${file_name} PRIVATE 
            OpenSSL::SSL 
            OpenSSL::Crypto
        )

    if(WIN32 AND file_name STREQUAL "test_consensus")
        # FHE kernels use large fixed-size temporaries.  Match the larger
        # pthread stack normally available on the Unix test hosts.
        target_link_options(${file_name} PRIVATE "-Wl,--stack,67108864")
    endif()
endforeach()
