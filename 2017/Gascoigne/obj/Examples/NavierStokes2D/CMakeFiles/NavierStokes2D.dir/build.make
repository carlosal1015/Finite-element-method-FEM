# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/carlosal1015/Gascoigne

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/carlosal1015/Gascoigne/obj

# Include any dependencies generated for this target.
include Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/depend.make

# Include the progress variables for this target.
include Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/progress.make

# Include the compile flags for this target's objects.
include Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/flags.make

Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/main.o: Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/flags.make
Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/main.o: ../Examples/NavierStokes2D/main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosal1015/Gascoigne/obj/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/main.o"
	cd /home/carlosal1015/Gascoigne/obj/Examples/NavierStokes2D && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NavierStokes2D.dir/main.o -c /home/carlosal1015/Gascoigne/Examples/NavierStokes2D/main.cc

Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/main.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NavierStokes2D.dir/main.i"
	cd /home/carlosal1015/Gascoigne/obj/Examples/NavierStokes2D && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosal1015/Gascoigne/Examples/NavierStokes2D/main.cc > CMakeFiles/NavierStokes2D.dir/main.i

Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/main.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NavierStokes2D.dir/main.s"
	cd /home/carlosal1015/Gascoigne/obj/Examples/NavierStokes2D && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosal1015/Gascoigne/Examples/NavierStokes2D/main.cc -o CMakeFiles/NavierStokes2D.dir/main.s

# Object files for target NavierStokes2D
NavierStokes2D_OBJECTS = \
"CMakeFiles/NavierStokes2D.dir/main.o"

# External object files for target NavierStokes2D
NavierStokes2D_EXTERNAL_OBJECTS =

Examples/NavierStokes2D/NavierStokes2D: Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/main.o
Examples/NavierStokes2D/NavierStokes2D: Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/build.make
Examples/NavierStokes2D/NavierStokes2D: ../lib/libGascoigneStd.so
Examples/NavierStokes2D/NavierStokes2D: Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/carlosal1015/Gascoigne/obj/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable NavierStokes2D"
	cd /home/carlosal1015/Gascoigne/obj/Examples/NavierStokes2D && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NavierStokes2D.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/build: Examples/NavierStokes2D/NavierStokes2D

.PHONY : Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/build

Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/clean:
	cd /home/carlosal1015/Gascoigne/obj/Examples/NavierStokes2D && $(CMAKE_COMMAND) -P CMakeFiles/NavierStokes2D.dir/cmake_clean.cmake
.PHONY : Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/clean

Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/depend:
	cd /home/carlosal1015/Gascoigne/obj && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/carlosal1015/Gascoigne /home/carlosal1015/Gascoigne/Examples/NavierStokes2D /home/carlosal1015/Gascoigne/obj /home/carlosal1015/Gascoigne/obj/Examples/NavierStokes2D /home/carlosal1015/Gascoigne/obj/Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Examples/NavierStokes2D/CMakeFiles/NavierStokes2D.dir/depend

