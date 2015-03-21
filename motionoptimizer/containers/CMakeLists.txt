set( SOURCE 
    ${CMAKE_CURRENT_SOURCE_DIR}/Trajectory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ChompGradient.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ConstraintFactory.cpp
   )

add_library(containers ${LIBRARY_TYPE}  ${SOURCE} )
target_link_libraries( containers mzcommon )

