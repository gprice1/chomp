
add_subdirectory( optimizer )
include_directories( nlopt )

add_library(chomp SHARED
    ChompOptimizerBase.cpp
    chomputil.cpp

    ConstraintFactory.cpp
    Constraint.cpp
    
    ChompOptimizer.cpp
    ChompGradient.cpp
    )

target_link_libraries( chomp mzcommon optimizer )

add_executable(testmatops testmatops.cpp)
target_link_libraries(testmatops chomp)

add_subdirectory(nlopt)
