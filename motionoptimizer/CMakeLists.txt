set(LINK_LIBS)
set( SOURCE )
set( HEADER )

add_subdirectory( utils )
include_directories( utils )

add_subdirectory( constraint )
include_directories( constraint )

add_subdirectory( containers )
include_directories( containers )

add_subdirectory( optimizer )
include_directories( optimizer )


set( SOURCE MotionOptimizer.cpp ${SOURCE} )
set( HEADER MotionOptimizer.h   ${HEADER} )

message(STATUS "Sources: ${SOURCE}")


add_library(motionoptimizer ${SOURCE})
target_link_libraries( motionoptimizer ${LINK_LIBS} )


