
# -*- mode: cmake; eval: (hl-line-mode 0);-*-

set(sources "main.cpp eigen.cpp interpolation.cpp")
set(linklibs "differential;gsl")

set(properties " ") 
build_and_install(${sources} "PROGRAM" ${linklibs} ${properties} "NOINSTALL")
