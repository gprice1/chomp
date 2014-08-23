add_subdirectory( utils )
include_directories( utils )

add_subdirectory( constraint )
include_directories( constraint )

add_subdirectory( optimizer )
include_directories( optimizer )

add_library(chomp SHARED Chomp.cpp )
target_link_libraries( chomp optimizer constraint utils mzcommon )


