
#use this to find the FindNLopt.cmake file.

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}")

find_package(NLopt)
include_directories(${NLOPT_INCLUDE_DIRS})
set(nlopt_libs ${NLOPT_LIBRARIES})


if ( NLOPT_FOUND )

    add_executable(nlopt_test nlopt_test.cpp)
    target_link_libraries( nlopt_test ${nlopt_libs} )

endif()




