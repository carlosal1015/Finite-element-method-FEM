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
include Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/depend.make

# Include the progress variables for this target.
include Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/progress.make

# Include the compile flags for this target's objects.
include Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/flags.make

Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/main1.o: Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/flags.make
Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/main1.o: ../Examples/Laplace2D/main1.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosal1015/Gascoigne/obj/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/main1.o"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Laplace2D && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Laplace2D-1.dir/main1.o -c /home/carlosal1015/Gascoigne/Examples/Laplace2D/main1.cc

Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/main1.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Laplace2D-1.dir/main1.i"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Laplace2D && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosal1015/Gascoigne/Examples/Laplace2D/main1.cc > CMakeFiles/Laplace2D-1.dir/main1.i

Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/main1.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Laplace2D-1.dir/main1.s"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Laplace2D && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosal1015/Gascoigne/Examples/Laplace2D/main1.cc -o CMakeFiles/Laplace2D-1.dir/main1.s

# Object files for target Laplace2D-1
Laplace2D__1_OBJECTS = \
"CMakeFiles/Laplace2D-1.dir/main1.o"

# External object files for target Laplace2D-1
Laplace2D__1_EXTERNAL_OBJECTS =

Examples/Laplace2D/Laplace2D-1: Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/main1.o
Examples/Laplace2D/Laplace2D-1: Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/build.make
Examples/Laplace2D/Laplace2D-1: ../lib/libGascoigneStd.so
Examples/Laplace2D/Laplace2D-1: Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/carlosal1015/Gascoigne/obj/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Laplace2D-1"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Laplace2D && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Laplace2D-1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/build: Examples/Laplace2D/Laplace2D-1

.PHONY : Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/build

Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/clean:
	cd /home/carlosal1015/Gascoigne/obj/Examples/Laplace2D && $(CMAKE_COMMAND) -P CMakeFiles/Laplace2D-1.dir/cmake_clean.cmake
.PHONY : Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/clean

Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/depend:
	cd /home/carlosal1015/Gascoigne/obj && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/carlosal1015/Gascoigne /home/carlosal1015/Gascoigne/Examples/Laplace2D /home/carlosal1015/Gascoigne/obj /home/carlosal1015/Gascoigne/obj/Examples/Laplace2D /home/carlosal1015/Gascoigne/obj/Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Examples/Laplace2D/CMakeFiles/Laplace2D-1.dir/depend

