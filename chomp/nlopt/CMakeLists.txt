

if ( NLOPT_FOUND )
    
    add_library( nlopt_chomp ChompNLopt.cpp )
    target_link_libraries( nlopt_chomp ${NLOPT_LIBRARIES} chomp ) 

    add_executable(nlopt_test nlopt_test.cpp)
    target_link_libraries( nlopt_test ${NLOPT_LIBRARIES} )

endif()

