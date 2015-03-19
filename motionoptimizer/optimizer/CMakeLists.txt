
set( LOCAL_SOURCE 
             
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizer.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompLocalOptimizer.cpp

     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizerBase.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/OptimizerBase.cpp
     
     ${CMAKE_CURRENT_SOURCE_DIR}/HMC.cpp
     ${SOURCE}
   )
 
set( LOCAL_HEADER 
    ${CMAKE_CURRENT_SOURCE_DIR}/OptimizerBase.h
             
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizerBase.h
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizer.h
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompLocalOptimizer.h

     ${CMAKE_CURRENT_SOURCE_DIR}/HMC.h
     ${HEADER}
   )


if ( NLOPT_FOUND )
    
    set( LOCAL_SOURCE 
         ${CMAKE_CURRENT_SOURCE_DIR}/NLOptimizer.cpp
         ${LOCAL_SOURCE} 
       )
    set( LOCAL_HEADER 
         ${CMAKE_CURRENT_SOURCE_DIR}/NLOptimizer.h
         ${LOCAL_HEADER} 
         )

    set( LINK_LIBS ${NLOPT_LIBRARIES}
         ${LINK_LIBS} 
         PARENT_SCOPE
       )

    add_executable(nlopt_test nlopt_test.cpp)
    target_link_libraries( nlopt_test ${NLOPT_LIBRARIES} )

endif( NLOPT_FOUND )

set( SOURCE 
     ${LOCAL_SOURCE}
     ${SOURCE}
     PARENT_SCOPE
     )
 set( HEADER 
     ${LOCAL_HEADER}
     ${HEADER}
     PARENT_SCOPE
     )
