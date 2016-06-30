# set minimum required cmake version
IF(APPLE)
  cmake_minimum_required (VERSION 2.8)
ELSE(APPLE)
  cmake_minimum_required (VERSION 2.6)
ENDIF(APPLE)


# set (BUILD_SHARED_LIBS TRUE)

# don't use configuration-specific subdirectories in Xcode
if (CMAKE_GENERATOR STREQUAL "Xcode")
   set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_SOURCE_DIR}/lib)
   set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}/lib)
endif (CMAKE_GENERATOR STREQUAL "Xcode")


# find Atom
find_file (ATOM_CONFIG atom-config PATHS ENV PATH NO_DEFAULT_PATH)
if (ATOM_CONFIG STREQUAL "ATOM_CONFIG-NOTFOUND")
   message (FATAL_ERROR "atom-config not found.  Make sure Atom is installed and \$(ATOM)/bin is in your path!")
else (ATOM_CONFIG STREQUAL "ATOM_CONFIG-NOTFOUND")
   message (STATUS "Found atom-config at ${ATOM_CONFIG}")
endif (ATOM_CONFIG STREQUAL "ATOM_CONFIG-NOTFOUND")

# get Atom flags
execute_process (COMMAND ${ATOM_CONFIG} --prefix
   OUTPUT_VARIABLE ATOM
)
execute_process (COMMAND ${ATOM_CONFIG} --cppflags
   OUTPUT_VARIABLE ATOM_CPPFLAGS
)
execute_process (COMMAND ${ATOM_CONFIG} --ldflags
   OUTPUT_VARIABLE ATOM_LDFLAGS
)
execute_process (COMMAND ${ATOM_CONFIG} --libs
   OUTPUT_VARIABLE ATOM_LIBS
)
string (STRIP ${ATOM} ATOM)
string (STRIP ${ATOM_CPPFLAGS} ATOM_CPPFLAGS)
string (STRIP ${ATOM_LDFLAGS} ATOM_LDFLAGS)
string (STRIP ${ATOM_LIBS} ATOM_LIBS)

set( CMAKE_CXX_FLAGS ${ATOM_CPPFLAGS} )

if (CMAKE_COMPILER_IS_GNUCXX)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -pedantic -Wall -Wextra -Wno-unused-parameter -Wno-long-long -std=c++11")
endif(CMAKE_COMPILER_IS_GNUCXX)

#setting fPIC
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

#setting correct rpaths
IF(APPLE)
  message(STATUS "Setting OSX rpath behavior.")
  set(CMAKE_MACOSX_RPATH 1)
ENDIF(APPLE)

# # find ROOT
# find_file (ROOT_CONFIG root-config PATHS ENV PATH NO_DEFAULT_PATH)
# if (ROOT_CONFIG STREQUAL "ROOT_CONFIG-NOTFOUND")
#    message (FATAL_ERROR "root-config not found.  Make sure Atom is installed and \$(ATOM)/bin is in your path!")
# else (ROOT_CONFIG STREQUAL "ROOT_CONFIG-NOTFOUND")
#    message (STATUS "Found root-config at ${ROOT_CONFIG}")
# endif (ROOT_CONFIG STREQUAL "ROOT_CONFIG-NOTFOUND")

# find_file (ROOTCINT rootcint PATHS ENV PATH NO_DEFAULT_PATH)
# if (ROOTCINT STREQUAL "ROOTCINT-NOTFOUND")
#    message (FATAL_ERROR "rootcint not found.  Make sure ROOT is installed and \$(ROOTSYS)/bin is in your path!")
# else (ROOTCINT STREQUAL "ROOTCINT-NOTFOUND")
#    message (STATUS "Found rootcint at ${ROOTCINT}")
# endif (ROOTCINT STREQUAL "ROOTCINT-NOTFOUND")


# # get ROOT flags
# execute_process (COMMAND ${ROOT_CONFIG} --prefix
#    OUTPUT_VARIABLE ROOTSYS
# )
# execute_process (COMMAND ${ROOT_CONFIG} --cflags
#    OUTPUT_VARIABLE ROOT_CPPFLAGS
# )
# execute_process (COMMAND ${ROOT_CONFIG} --ldflags
#    OUTPUT_VARIABLE ROOT_LDFLAGS
# )
# execute_process (COMMAND ${ROOT_CONFIG} --libs
#    OUTPUT_VARIABLE ROOT_LIBS
# )
# string (STRIP ${ROOTSYS} ROOTSYS)
# string (STRIP ${ROOT_CPPFLAGS} ROOT_CPPFLAGS)
# string (STRIP ${ROOT_LDFLAGS} ROOT_LDFLAGS)
# string (STRIP ${ROOT_LIBS} ROOT_LIBS)


# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CPPFLAGS}")

set (Atom_Libraries AtomBase AtomCore AtomProjections AtomSimulation DetectorYaml AtomTools RivetTools /usr/local/lib/libHepMC.dylib -Wl,-rpath,/usr/local/Cellar/fastjet/3.1.3/lib -lm -L/usr/local/Cellar/fastjet/3.1.3/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone /usr/local/Cellar/yoda/1.5.9/lib/libYODA.dylib)
function (new_ana name)
  add_library(Atom_${name} SHARED ${name}.cc)
  file(READ ${name}.cc Atom_${name}_source)
  string(REGEX REPLACE " |\n|\t" ""  Atom_${name}_repl ${Atom_${name}_source})
  string(MD5 Atom_${name}_md5 ${Atom_${name}_repl})
  file(APPEND "${PROJECT_SOURCE_DIR}/Analyses_hash.txt" "${name}: ${Atom_${name}_md5}\n")
  if(${ARGC} GREATER 1) 
	  target_link_libraries(Atom_${name} ${Atom_Libraries} ${ARGV1})
    get_property(new_link_flags TARGET Atom_${name} PROPERTY LINK_FLAGS)
    set(new_link_flags "${new_link_flags} ${ATOM_LDFLAGS}")
    SET_TARGET_PROPERTIES(Atom_${name} PROPERTIES LINK_FLAGS ${new_link_flags})
  else()
	  target_link_libraries(Atom_${name} ${Atom_Libraries})
    get_property(new_link_flags TARGET Atom_${name} PROPERTY LINK_FLAGS)
    set(new_link_flags "${new_link_flags} ${ATOM_LDFLAGS}")
    SET_TARGET_PROPERTIES(Atom_${name} PROPERTIES LINK_FLAGS ${new_link_flags})
  endif(${ARGC} GREATER 1)
  file(READ ${ATOM}/bin/atom-addanalysispath.csh AddAnalysisPathCsh_source)
  string(REPLACE "$1" "${PROJECT_SOURCE_DIR}" AddAnalysisPathCsh_repl ${AddAnalysisPathCsh_source})
  file(WRITE "${PROJECT_SOURCE_DIR}/addAnalysisPath.csh" ${AddAnalysisPathCsh_repl})
  file(READ ${ATOM}/bin/atom-addanalysispath.sh AddAnalysisPathSh_source)
  string(REPLACE "$1" "${PROJECT_SOURCE_DIR}" AddAnalysisPathSh_repl ${AddAnalysisPathSh_source})
  file(WRITE "${PROJECT_SOURCE_DIR}/addAnalysisPath.sh" ${AddAnalysisPathSh_repl})
endfunction (new_ana)

file(GLOB CustomAnalyses "${PROJECT_SOURCE_DIR}/*.cc") 
foreach(anasourcefile ${CustomAnalyses})
  STRING(REGEX REPLACE "(.*)/(.*).cc" "\\2" ananame "${anasourcefile}")
  MESSAGE(STATUS "Process Analysis ${ananame}")
  new_ana(${ananame})
endforeach()
