# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_SOURCE_DIR = /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build

# Include any dependencies generated for this target.
include test/CMakeFiles/kruskal_test.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/kruskal_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/kruskal_test.dir/flags.make

test/CMakeFiles/kruskal_test.dir/kruskal_test.cc.o: test/CMakeFiles/kruskal_test.dir/flags.make
test/CMakeFiles/kruskal_test.dir/kruskal_test.cc.o: ../test/kruskal_test.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/kruskal_test.dir/kruskal_test.cc.o"
	cd /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/kruskal_test.dir/kruskal_test.cc.o -c /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/test/kruskal_test.cc

test/CMakeFiles/kruskal_test.dir/kruskal_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/kruskal_test.dir/kruskal_test.cc.i"
	cd /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/test/kruskal_test.cc > CMakeFiles/kruskal_test.dir/kruskal_test.cc.i

test/CMakeFiles/kruskal_test.dir/kruskal_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/kruskal_test.dir/kruskal_test.cc.s"
	cd /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/test/kruskal_test.cc -o CMakeFiles/kruskal_test.dir/kruskal_test.cc.s

# Object files for target kruskal_test
kruskal_test_OBJECTS = \
"CMakeFiles/kruskal_test.dir/kruskal_test.cc.o"

# External object files for target kruskal_test
kruskal_test_EXTERNAL_OBJECTS =

test/kruskal_test: test/CMakeFiles/kruskal_test.dir/kruskal_test.cc.o
test/kruskal_test: test/CMakeFiles/kruskal_test.dir/build.make
test/kruskal_test: lemon/libemon.a
test/kruskal_test: test/CMakeFiles/kruskal_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable kruskal_test"
	cd /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/kruskal_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/kruskal_test.dir/build: test/kruskal_test

.PHONY : test/CMakeFiles/kruskal_test.dir/build

test/CMakeFiles/kruskal_test.dir/clean:
	cd /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -P CMakeFiles/kruskal_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/kruskal_test.dir/clean

test/CMakeFiles/kruskal_test.dir/depend:
	cd /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1 /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/test /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/test /home/daniel/Dropbox/kth/kth-dd2467-individual-project/clustering/lemon-1.3.1/build/test/CMakeFiles/kruskal_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/kruskal_test.dir/depend

