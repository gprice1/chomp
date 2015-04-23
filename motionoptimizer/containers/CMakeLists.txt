set( SOURCE 
    ${CMAKE_CURRENT_SOURCE_DIR}/Trajectory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Gradient.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ConstraintFactory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ProblemDescription.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Metric.cpp
   )

add_library(containers ${LIBRARY_TYPE} ${SOURCE} )
target_link_libraries( containers mzcommon )

add_executable(testmetric testmetric.cpp )
target_link_libraries( testmetric containers ) 

