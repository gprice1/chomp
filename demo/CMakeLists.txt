
set(demo_libs mzcommon motionoptimizer ${PNG_LIBRARY})

if (${CAIRO_FOUND})
    set(demo_libs ${demo_libs} ${CAIRO_LDFLAGS})
endif()

add_executable(circle_demo circle_demo.cpp)
target_link_libraries(circle_demo ${demo_libs})


#add_executable(circle_goalset_demo circle_goalset_demo.cpp)
#target_link_libraries(circle_goalset_demo ${demo_libs})

add_executable(map2d_demo map2d_demo.cpp Map2D.cpp)
target_link_libraries(map2d_demo ${demo_libs})

add_executable(map2d_eval map2d_eval.cpp Map2D.cpp)
target_link_libraries(map2d_eval ${demo_libs})

