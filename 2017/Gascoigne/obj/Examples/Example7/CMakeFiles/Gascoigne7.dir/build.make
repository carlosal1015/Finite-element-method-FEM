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
include Examples/Example7/CMakeFiles/Gascoigne7.dir/depend.make

# Include the progress variables for this target.
include Examples/Example7/CMakeFiles/Gascoigne7.dir/progress.make

# Include the compile flags for this target's objects.
include Examples/Example7/CMakeFiles/Gascoigne7.dir/flags.make

Examples/Example7/CMakeFiles/Gascoigne7.dir/dwralgorithm.o: Examples/Example7/CMakeFiles/Gascoigne7.dir/flags.make
Examples/Example7/CMakeFiles/Gascoigne7.dir/dwralgorithm.o: ../Examples/Example7/dwralgorithm.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosal1015/Gascoigne/obj/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Examples/Example7/CMakeFiles/Gascoigne7.dir/dwralgorithm.o"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Example7 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Gascoigne7.dir/dwralgorithm.o -c /home/carlosal1015/Gascoigne/Examples/Example7/dwralgorithm.cc

Examples/Example7/CMakeFiles/Gascoigne7.dir/dwralgorithm.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Gascoigne7.dir/dwralgorithm.i"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Example7 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosal1015/Gascoigne/Examples/Example7/dwralgorithm.cc > CMakeFiles/Gascoigne7.dir/dwralgorithm.i

Examples/Example7/CMakeFiles/Gascoigne7.dir/dwralgorithm.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Gascoigne7.dir/dwralgorithm.s"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Example7 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosal1015/Gascoigne/Examples/Example7/dwralgorithm.cc -o CMakeFiles/Gascoigne7.dir/dwralgorithm.s

Examples/Example7/CMakeFiles/Gascoigne7.dir/main.o: Examples/Example7/CMakeFiles/Gascoigne7.dir/flags.make
Examples/Example7/CMakeFiles/Gascoigne7.dir/main.o: ../Examples/Example7/main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosal1015/Gascoigne/obj/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Examples/Example7/CMakeFiles/Gascoigne7.dir/main.o"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Example7 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Gascoigne7.dir/main.o -c /home/carlosal1015/Gascoigne/Examples/Example7/main.cc

Examples/Example7/CMakeFiles/Gascoigne7.dir/main.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Gascoigne7.dir/main.i"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Example7 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosal1015/Gascoigne/Examples/Example7/main.cc > CMakeFiles/Gascoigne7.dir/main.i

Examples/Example7/CMakeFiles/Gascoigne7.dir/main.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Gascoigne7.dir/main.s"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Example7 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosal1015/Gascoigne/Examples/Example7/main.cc -o CMakeFiles/Gascoigne7.dir/main.s

# Object files for target Gascoigne7
Gascoigne7_OBJECTS = \
"CMakeFiles/Gascoigne7.dir/dwralgorithm.o" \
"CMakeFiles/Gascoigne7.dir/main.o"

# External object files for target Gascoigne7
Gascoigne7_EXTERNAL_OBJECTS =

Examples/Example7/Gascoigne7: Examples/Example7/CMakeFiles/Gascoigne7.dir/dwralgorithm.o
Examples/Example7/Gascoigne7: Examples/Example7/CMakeFiles/Gascoigne7.dir/main.o
Examples/Example7/Gascoigne7: Examples/Example7/CMakeFiles/Gascoigne7.dir/build.make
Examples/Example7/Gascoigne7: ../lib/libGascoigneStd.so
Examples/Example7/Gascoigne7: Examples/Example7/CMakeFiles/Gascoigne7.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/carlosal1015/Gascoigne/obj/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable Gascoigne7"
	cd /home/carlosal1015/Gascoigne/obj/Examples/Example7 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Gascoigne7.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Examples/Example7/CMakeFiles/Gascoigne7.dir/build: Examples/Example7/Gascoigne7

.PHONY : Examples/Example7/CMakeFiles/Gascoigne7.dir/build

Examples/Example7/CMakeFiles/Gascoigne7.dir/clean:
	cd /home/carlosal1015/Gascoigne/obj/Examples/Example7 && $(CMAKE_COMMAND) -P CMakeFiles/Gascoigne7.dir/cmake_clean.cmake
.PHONY : Examples/Example7/CMakeFiles/Gascoigne7.dir/clean

Examples/Example7/CMakeFiles/Gascoigne7.dir/depend:
	cd /home/carlosal1015/Gascoigne/obj && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/carlosal1015/Gascoigne /home/carlosal1015/Gascoigne/Examples/Example7 /home/carlosal1015/Gascoigne/obj /home/carlosal1015/Gascoigne/obj/Examples/Example7 /home/carlosal1015/Gascoigne/obj/Examples/Example7/CMakeFiles/Gascoigne7.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Examples/Example7/CMakeFiles/Gascoigne7.dir/depend

