set( LINK_LIBS )

set( SOURCE 
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizer.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompLocalOptimizer.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizerBase.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/OptimizerBase.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/HMC.cpp
   )
 

if ( NLOPT_FOUND )
    
    set( SOURCE 
         ${CMAKE_CURRENT_SOURCE_DIR}/NLOptimizer.cpp
         ${SOURCE}
       )

    set( LINK_LIBS ${NLOPT_LIBRARIES} )

    add_executable(nlopt_test nlopt_test.cpp)
    target_link_libraries( nlopt_test ${NLOPT_LIBRARIES} )

endif( NLOPT_FOUND )

add_library( optimizer SHARED ${SOURCE} )
target_link_libraries( optimizer ${LINK_LIBS} containers)
