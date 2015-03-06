
set( SOURCE ${SOURCE}
     ${CMAKE_CURRENT_SOURCE_DIR}/OptimizerBase.cpp
             
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizerBase.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizer.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompLocalOptimizer.cpp

     ${CMAKE_CURRENT_SOURCE_DIR}/HMC.cpp
     PARENT_SCOPE
   )
 
set( HEADER ${HEADER}
     ${CMAKE_CURRENT_SOURCE_DIR}/OptimizerBase.h
             
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizerBase.h
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizer.h
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompLocalOptimizer.h

     ${CMAKE_CURRENT_SOURCE_DIR}/HMC.h
     PARENT_SCOPE
   )


if ( NLOPT_FOUND )
    
    set( SOURCE ${SOURCE} 
         ${CMAKE_CURRENT_SOURCE_DIR}/NLOptimizer.cpp
         PARENT_SCOPE
       )
    set( LINK_LIBS ${NLOPT_LIBRARIES}
         ${LINK_LIBS} 
         PARENT_SCOPE
       )

    add_executable(nlopt_test nlopt_test.cpp)
    target_link_libraries( nlopt_test ${NLOPT_LIBRARIES} )

endif( NLOPT_FOUND )
