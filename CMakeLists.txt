# -*- mode: cmake; eval: (hl-line-mode 0);-*-
cmake_minimum_required(VERSION 3.0.0)
project(int VERSION 1.0 LANGUAGES CXX)

set(libs resolv)
add_compile_options(-pipe -O2 -Wall -W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fno-common --std=gnu++17)

set(progname ${PROJECT_NAME})
add_definitions(-Dcmake_MAJOR=${PROJECT_VERSION_MAJOR} -Dcmake_MINOR=${PROJECT_VERSION_MINOR} -Dcmake_NAME="${progname}" -D___debug="1")

include("default_buildinstall.cmake")

foreach(loopvar IN LISTS libs)
  message("entering ${loopvar}")
  include_directories(${loopvar})
  add_subdirectory(${loopvar} ${CMAKE_CURRENT_BINARY_DIR}/${loopvar})
  message("${INCLUDE_DIRECTORIES}")
  message("${SYSTEM_INCLUDE_DIRECTORIES}")
  message("${CMAKE_INCLUDE_DIRECTORIES}")
  message("leaving ${loopvar}")
endforeach(loopvar)

install(FILES LICENSE DESTINATION /usr/local/share/${progname} PERMISSIONS WORLD_READ)
