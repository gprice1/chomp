

add_library(chomp SHARED
    ChompOptimizerBase.cpp
    chomputil.cpp

    ConstraintFactory.cpp
    Constraint.cpp
    
    ChompOptimizer.cpp
    ChompGradient.cpp
    
    HMC.cpp
    )

target_link_libraries( chomp mzcommon )

add_executable(testmatops testmatops.cpp)
target_link_libraries(testmatops chomp)

add_subdirectory(nlopt)
include_directories( nlopt )
