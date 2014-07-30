

add_library(chomp SHARED
    chomputil.cpp
    
    ConstraintFactory.cpp
    Constraint.cpp
    
    Chomp.cpp
    ChompGradient.cpp
    
    HMC.cpp
    )



target_link_libraries( chomp mzcommon )

add_executable(testmatops testmatops.cpp)
target_link_libraries(testmatops chomp)

add_executable(circle_demo circle_demo.cpp)
add_executable(circle_goalset_demo circle_goalset_demo.cpp)

add_executable(map2d_demo map2d_demo.cpp Map2D.cpp)

set(demo_libs mzcommon chomp ${PNG_LIBRARY})

if (${CAIRO_FOUND})
set(demo_libs ${demo_libs} ${CAIRO_LDFLAGS})
endif()

target_link_libraries(circle_goalset_demo ${demo_libs})
target_link_libraries(circle_demo ${demo_libs})
target_link_libraries(map2d_demo ${demo_libs})


