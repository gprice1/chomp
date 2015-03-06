set( SOURCE MotionOptimizer.cpp )
set( HEADER MotionOptimizer.h  )
set(LINK_LIBS)

add_subdirectory( utils )
include_directories( utils )

add_subdirectory( constraint )
include_directories( constraint )

add_subdirectory( containers )
include_directories( containers )

add_subdirectory( optimizer )
include_directories( optimizer )

add_library(motion_optimizer ${SOURCE})
target_link_libraries( motion_optimizer ${LINK_LIBS} )


