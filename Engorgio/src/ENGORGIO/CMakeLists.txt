set(ENGORGIO_SOURCE_FILES ${ENGORGIO_SOURCE_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/comp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/quant.cpp
    ${CMAKE_CURRENT_LIST_DIR}/orderapp.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ordergen.cpp
    ${CMAKE_CURRENT_LIST_DIR}/utils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/sort.cpp
)

# Add header files for installation
install(
    FILES
        ${CMAKE_CURRENT_LIST_DIR}/comp.h
        ${CMAKE_CURRENT_LIST_DIR}/quant.h 
        ${CMAKE_CURRENT_LIST_DIR}/orderapp.h
        ${CMAKE_CURRENT_LIST_DIR}/ordergen.h 
        ${CMAKE_CURRENT_LIST_DIR}/utils.h       
        ${CMAKE_CURRENT_LIST_DIR}/sort.h
    DESTINATION
    ${ENGORGIO_INCLUDES_INSTALL_DIR}/ENGORGIO
)

set(ENGORGIO_SOURCE_FILES ${ENGORGIO_SOURCE_FILES} PARENT_SCOPE)