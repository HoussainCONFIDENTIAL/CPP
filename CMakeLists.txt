cmake_minimum_required(VERSION 2.6) # Compatibility with CMake
project(create_raster) # Project name: "create_raster" 

set(CMAKE_BUILD_TYPE debug) # Compilation in debug mode
set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++11") # Options

file( # List of files to compile
  GLOB_RECURSE # Recursive search  
  source_files # Files listed in a variable source_files  
  src/*.cpp # Form of files to search  
)

add_executable( # Creation of the executable binary for the project 
  ${PROJECT_NAME} # Contains the name of the binary  
  ${source_files} # List of files compiled for linking
        )

target_link_libraries(${PROJECT_NAME} proj)

# Inorder to setup the "create_raster" project, link the PROJ library to it, and add an excutable named just like the project, I wrote this CMakeLists.txt.
