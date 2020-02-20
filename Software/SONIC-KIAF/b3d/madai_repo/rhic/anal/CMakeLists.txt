project( analyze )
cmake_minimum_required( VERSION 2.8.0 )

# These options only need to be set if building this project by itself
if ( NOT rhic_BINARY_DIR )
  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )

  ## Need coral library
  find_package( coral REQUIRED )
endif()

set( analyze_INCLUDE_DIRS
  ${analyze_SOURCE_DIR}/src
  ${coral_INCLUDE_DIRS}
)
set( analyze_LIBRARY_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} )
set( analyze_LIBRARIES analyze ${coral_LIBRARIES} )

add_subdirectory( src )

## Configure the package file
configure_file( analyze-config.cmake.in ${PROJECT_BINARY_DIR}/analyze-config.cmake @ONLY )
