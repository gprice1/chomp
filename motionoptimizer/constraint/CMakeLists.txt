set( SOURCE 
    ${CMAKE_CURRENT_SOURCE_DIR}/ConstraintFactory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Constraint.cpp
   )

add_library( constraint SHARED ${SOURCE} )
target_link_libraries( constraint mzcommon )
