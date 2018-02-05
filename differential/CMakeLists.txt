# -*- mode: cmake; eval: (hl-line-mode 0);-*-
set(sources "differential.cpp")
set(linklibs " ")
set(properties " ") 
build_and_install(${sources} "SHARED" ${linklibs} ${properties} "NOINSTALL")
