# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.5.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.5.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054

# Include any dependencies generated for this target.
include CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/flags.make

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o: CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/flags.make
CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o: ATLAS_CONF_2016_054.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o -c /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054/ATLAS_CONF_2016_054.cc

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054/ATLAS_CONF_2016_054.cc > CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.i

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054/ATLAS_CONF_2016_054.cc -o CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.s

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o.requires:

.PHONY : CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o.requires

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o.provides: CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o.requires
	$(MAKE) -f CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/build.make CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o.provides.build
.PHONY : CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o.provides

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o.provides.build: CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o


# Object files for target Atom_ATLAS_CONF_2016_054
Atom_ATLAS_CONF_2016_054_OBJECTS = \
"CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o"

# External object files for target Atom_ATLAS_CONF_2016_054
Atom_ATLAS_CONF_2016_054_EXTERNAL_OBJECTS =

libAtom_ATLAS_CONF_2016_054.dylib: CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o
libAtom_ATLAS_CONF_2016_054.dylib: CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/build.make
libAtom_ATLAS_CONF_2016_054.dylib: /usr/local/lib/libHepMC.dylib
libAtom_ATLAS_CONF_2016_054.dylib: /usr/local/Cellar/yoda/1.5.9/lib/libYODA.dylib
libAtom_ATLAS_CONF_2016_054.dylib: CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libAtom_ATLAS_CONF_2016_054.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/build: libAtom_ATLAS_CONF_2016_054.dylib

.PHONY : CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/build

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/requires: CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/ATLAS_CONF_2016_054.cc.o.requires

.PHONY : CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/requires

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/clean

CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/depend:
	cd /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054 /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054 /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054 /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054 /Users/spliew/Desktop/tools/Atom_src/Atom-work/analyses/ATLAS-CONF-2016-054/CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Atom_ATLAS_CONF_2016_054.dir/depend

