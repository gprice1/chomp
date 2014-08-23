add_library(utils
            ChompGradient.cpp
            Trajectory.cpp
            )


add_executable(testmatops testmatops.cpp chomputil.h)

