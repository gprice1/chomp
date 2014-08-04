

add_library(chomp SHARED
    ChompOptimizerBase.cpp
    chomputil.cpp

    ConstraintFactory.cpp
    Constraint.cpp
    
    Chomp.cpp
    ChompGradient.cpp
    
    HMC.cpp
    )



target_link_libraries( chomp mzcommon )
