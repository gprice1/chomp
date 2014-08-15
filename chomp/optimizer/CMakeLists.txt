
set( sources 
     OptimizerBase.cpp
             
     ChompOptimizerBase.cpp
     ChompOptimizer.cpp
     ChompLocalOptimizer.cpp
     HMC.cpp
     )

set( link_libs )

if ( NLOPT_FOUND )
    
    set( sources ${sources} NLOptimizer.cpp )
    set( link_libs ${NLOPT_LIBRARIES} ${link_libs} )

    add_executable(nlopt_test nlopt_test.cpp)
    target_link_libraries( nlopt_test ${NLOPT_LIBRARIES} )

endif( NLOPT_FOUND )

add_library( optimizer SHARED ${sources} ) 
target_link_libraries( optimizer ${link_libs} )
