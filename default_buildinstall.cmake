# -*- mode: cmake; eval: (hl-line-mode 0);-*-
function(build_and_install sources type linklibs properties install)
  string(REPLACE "${CMAKE_SOURCE_DIR}/" "" progname_tmp ${CMAKE_CURRENT_SOURCE_DIR})
  string(REPLACE "/" "_" progname ${progname_tmp})
  set(progname ${progname} PARENT_SCOPE)

  if(${type} STREQUAL "PROGRAM")
    message("program")
    add_executable(${progname} ${sources})
  elseif(${type} STREQUAL "SHARED")
    message("lib")
    add_library(${progname} SHARED ${sources})
  else()
    message("Target ${progname} won't build anything")
    return()
  endif()
  foreach(loopvar IN LISTS properties)    
    set_property(TARGET ${progname} PROPERTY ${loopvar})
  endforeach(loopvar)
  target_link_libraries(${progname} ${linklibs})
  if(${install} MATCHES "NO.*")
    message("Target ${progname} won't install anything")
  else()
    if(${type} STREQUAL "PROGRAM")
      install(TARGETS ${progname} RUNTIME DESTINATION "/usr/local/bin/")
    elseif(${type} STREQUAL "SHARED")
      install(TARGETS ${progname} LIBRARY DESTINATION "/usr/local/lib/")
    endif()
  endif()
endfunction()
