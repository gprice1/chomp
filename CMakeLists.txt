# CMakeLists.txt for hubomz
cmake_minimum_required(VERSION 2.6)

set(CMAKE_BUILD_TYPE Debug)

set(BUILD_DEMOS TRUE)
set(BUILD_TESTS FALSE)

project(CHOMP_MULTIGRID)

include(mzsettings.cmake)

set(BASE_LIBS
  ${OPENGL_LIBRARY}
)

set(GUI_LIBS
  ${BASE_LIBS}
  ${GLUT_LIBRARY}
  ${PNG_LIBRARY}
)

add_subdirectory(mzcommon)
include_directories(mzcommon)

add_subdirectory(motionoptimizer)
include_directories( motionoptimizer )

if( BUILD_DEMOS )
    add_subdirectory( demo )
endif( BUILD_DEMOS )

link_directories ( ${CMAKE_CURRENT_SOURCE_DIR}/mzcommon
                   ${CMAKE_CURRENT_SOURCE_DIR}/motionoptimizer
                   )
