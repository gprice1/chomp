set( LINK_LIBS )

set( SOURCE 
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizer.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompLocalOptimizer.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/ChompOptimizerBase.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/TestOptimizer.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/OptimizerBase.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/HMC.cpp
   )
 

if ( NLOPT_FOUND AND USE_NLOPT )
    
    set( SOURCE 
         ${CMAKE_CURRENT_SOURCE_DIR}/NLOptimizer.cpp
         ${SOURCE}
       )

    set( LINK_LIBS ${NLOPT_LIBRARIES} )
    
    if( BUILD_TESTS )
        add_executable(nlopt_test nlopt_test.cpp)
        target_link_libraries( nlopt_test ${NLOPT_LIBRARIES} )
    endif( BUILD_TESTS )

endif()

add_library( optimizer ${LIBRARY_TYPE}  ${SOURCE} )
target_link_libraries( optimizer ${LINK_LIBS} containers)
