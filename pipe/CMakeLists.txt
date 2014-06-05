add_library( pipe Pipe.cpp)

add_executable( testpipe testpipe.cpp)
target_link_libraries( testpipe pipe chomp )


