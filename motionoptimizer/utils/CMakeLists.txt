
set( HEADER ${HEADER}
    ${CMAKE_CURRENT_SOURCE_DIR}/function_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/class_utils.h
    PARENT_SCOPE
   )

add_executable(testmatops testmatops.cpp function_utils.h class_utils.h )
target_link_libraries( testmatops mzcommon ) 

