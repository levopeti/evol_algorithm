# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /home/biot/clion-2018.1.2/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/biot/clion-2018.1.2/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/biot/projects/szakdolgozat/Evolutionary_algorithm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/untitled.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/untitled.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/untitled.dir/flags.make

CMakeFiles/untitled.dir/main.cpp.o: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/untitled.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/untitled.dir/main.cpp.o -c /home/biot/projects/szakdolgozat/Evolutionary_algorithm/main.cpp

CMakeFiles/untitled.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/biot/projects/szakdolgozat/Evolutionary_algorithm/main.cpp > CMakeFiles/untitled.dir/main.cpp.i

CMakeFiles/untitled.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/biot/projects/szakdolgozat/Evolutionary_algorithm/main.cpp -o CMakeFiles/untitled.dir/main.cpp.s

CMakeFiles/untitled.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/untitled.dir/main.cpp.o.requires

CMakeFiles/untitled.dir/main.cpp.o.provides: CMakeFiles/untitled.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/untitled.dir/build.make CMakeFiles/untitled.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/untitled.dir/main.cpp.o.provides

CMakeFiles/untitled.dir/main.cpp.o.provides.build: CMakeFiles/untitled.dir/main.cpp.o


CMakeFiles/untitled.dir/Gene.cpp.o: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Gene.cpp.o: ../Gene.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/untitled.dir/Gene.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/untitled.dir/Gene.cpp.o -c /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Gene.cpp

CMakeFiles/untitled.dir/Gene.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Gene.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Gene.cpp > CMakeFiles/untitled.dir/Gene.cpp.i

CMakeFiles/untitled.dir/Gene.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Gene.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Gene.cpp -o CMakeFiles/untitled.dir/Gene.cpp.s

CMakeFiles/untitled.dir/Gene.cpp.o.requires:

.PHONY : CMakeFiles/untitled.dir/Gene.cpp.o.requires

CMakeFiles/untitled.dir/Gene.cpp.o.provides: CMakeFiles/untitled.dir/Gene.cpp.o.requires
	$(MAKE) -f CMakeFiles/untitled.dir/build.make CMakeFiles/untitled.dir/Gene.cpp.o.provides.build
.PHONY : CMakeFiles/untitled.dir/Gene.cpp.o.provides

CMakeFiles/untitled.dir/Gene.cpp.o.provides.build: CMakeFiles/untitled.dir/Gene.cpp.o


CMakeFiles/untitled.dir/Swarm.cpp.o: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Swarm.cpp.o: ../Swarm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/untitled.dir/Swarm.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/untitled.dir/Swarm.cpp.o -c /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Swarm.cpp

CMakeFiles/untitled.dir/Swarm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Swarm.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Swarm.cpp > CMakeFiles/untitled.dir/Swarm.cpp.i

CMakeFiles/untitled.dir/Swarm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Swarm.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Swarm.cpp -o CMakeFiles/untitled.dir/Swarm.cpp.s

CMakeFiles/untitled.dir/Swarm.cpp.o.requires:

.PHONY : CMakeFiles/untitled.dir/Swarm.cpp.o.requires

CMakeFiles/untitled.dir/Swarm.cpp.o.provides: CMakeFiles/untitled.dir/Swarm.cpp.o.requires
	$(MAKE) -f CMakeFiles/untitled.dir/build.make CMakeFiles/untitled.dir/Swarm.cpp.o.provides.build
.PHONY : CMakeFiles/untitled.dir/Swarm.cpp.o.provides

CMakeFiles/untitled.dir/Swarm.cpp.o.provides.build: CMakeFiles/untitled.dir/Swarm.cpp.o


CMakeFiles/untitled.dir/Iteration.cpp.o: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Iteration.cpp.o: ../Iteration.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/untitled.dir/Iteration.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/untitled.dir/Iteration.cpp.o -c /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Iteration.cpp

CMakeFiles/untitled.dir/Iteration.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Iteration.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Iteration.cpp > CMakeFiles/untitled.dir/Iteration.cpp.i

CMakeFiles/untitled.dir/Iteration.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Iteration.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Iteration.cpp -o CMakeFiles/untitled.dir/Iteration.cpp.s

CMakeFiles/untitled.dir/Iteration.cpp.o.requires:

.PHONY : CMakeFiles/untitled.dir/Iteration.cpp.o.requires

CMakeFiles/untitled.dir/Iteration.cpp.o.provides: CMakeFiles/untitled.dir/Iteration.cpp.o.requires
	$(MAKE) -f CMakeFiles/untitled.dir/build.make CMakeFiles/untitled.dir/Iteration.cpp.o.provides.build
.PHONY : CMakeFiles/untitled.dir/Iteration.cpp.o.provides

CMakeFiles/untitled.dir/Iteration.cpp.o.provides.build: CMakeFiles/untitled.dir/Iteration.cpp.o


CMakeFiles/untitled.dir/Population.cpp.o: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Population.cpp.o: ../Population.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/untitled.dir/Population.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/untitled.dir/Population.cpp.o -c /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Population.cpp

CMakeFiles/untitled.dir/Population.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Population.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Population.cpp > CMakeFiles/untitled.dir/Population.cpp.i

CMakeFiles/untitled.dir/Population.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Population.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Population.cpp -o CMakeFiles/untitled.dir/Population.cpp.s

CMakeFiles/untitled.dir/Population.cpp.o.requires:

.PHONY : CMakeFiles/untitled.dir/Population.cpp.o.requires

CMakeFiles/untitled.dir/Population.cpp.o.provides: CMakeFiles/untitled.dir/Population.cpp.o.requires
	$(MAKE) -f CMakeFiles/untitled.dir/build.make CMakeFiles/untitled.dir/Population.cpp.o.provides.build
.PHONY : CMakeFiles/untitled.dir/Population.cpp.o.provides

CMakeFiles/untitled.dir/Population.cpp.o.provides.build: CMakeFiles/untitled.dir/Population.cpp.o


CMakeFiles/untitled.dir/Chromosome.cpp.o: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/Chromosome.cpp.o: ../Chromosome.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/untitled.dir/Chromosome.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/untitled.dir/Chromosome.cpp.o -c /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Chromosome.cpp

CMakeFiles/untitled.dir/Chromosome.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/Chromosome.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Chromosome.cpp > CMakeFiles/untitled.dir/Chromosome.cpp.i

CMakeFiles/untitled.dir/Chromosome.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/Chromosome.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/biot/projects/szakdolgozat/Evolutionary_algorithm/Chromosome.cpp -o CMakeFiles/untitled.dir/Chromosome.cpp.s

CMakeFiles/untitled.dir/Chromosome.cpp.o.requires:

.PHONY : CMakeFiles/untitled.dir/Chromosome.cpp.o.requires

CMakeFiles/untitled.dir/Chromosome.cpp.o.provides: CMakeFiles/untitled.dir/Chromosome.cpp.o.requires
	$(MAKE) -f CMakeFiles/untitled.dir/build.make CMakeFiles/untitled.dir/Chromosome.cpp.o.provides.build
.PHONY : CMakeFiles/untitled.dir/Chromosome.cpp.o.provides

CMakeFiles/untitled.dir/Chromosome.cpp.o.provides.build: CMakeFiles/untitled.dir/Chromosome.cpp.o


CMakeFiles/untitled.dir/python_test.cpp.o: CMakeFiles/untitled.dir/flags.make
CMakeFiles/untitled.dir/python_test.cpp.o: ../python_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/untitled.dir/python_test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/untitled.dir/python_test.cpp.o -c /home/biot/projects/szakdolgozat/Evolutionary_algorithm/python_test.cpp

CMakeFiles/untitled.dir/python_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/untitled.dir/python_test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/biot/projects/szakdolgozat/Evolutionary_algorithm/python_test.cpp > CMakeFiles/untitled.dir/python_test.cpp.i

CMakeFiles/untitled.dir/python_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/untitled.dir/python_test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/biot/projects/szakdolgozat/Evolutionary_algorithm/python_test.cpp -o CMakeFiles/untitled.dir/python_test.cpp.s

CMakeFiles/untitled.dir/python_test.cpp.o.requires:

.PHONY : CMakeFiles/untitled.dir/python_test.cpp.o.requires

CMakeFiles/untitled.dir/python_test.cpp.o.provides: CMakeFiles/untitled.dir/python_test.cpp.o.requires
	$(MAKE) -f CMakeFiles/untitled.dir/build.make CMakeFiles/untitled.dir/python_test.cpp.o.provides.build
.PHONY : CMakeFiles/untitled.dir/python_test.cpp.o.provides

CMakeFiles/untitled.dir/python_test.cpp.o.provides.build: CMakeFiles/untitled.dir/python_test.cpp.o


# Object files for target untitled
untitled_OBJECTS = \
"CMakeFiles/untitled.dir/main.cpp.o" \
"CMakeFiles/untitled.dir/Gene.cpp.o" \
"CMakeFiles/untitled.dir/Swarm.cpp.o" \
"CMakeFiles/untitled.dir/Iteration.cpp.o" \
"CMakeFiles/untitled.dir/Population.cpp.o" \
"CMakeFiles/untitled.dir/Chromosome.cpp.o" \
"CMakeFiles/untitled.dir/python_test.cpp.o"

# External object files for target untitled
untitled_EXTERNAL_OBJECTS =

untitled: CMakeFiles/untitled.dir/main.cpp.o
untitled: CMakeFiles/untitled.dir/Gene.cpp.o
untitled: CMakeFiles/untitled.dir/Swarm.cpp.o
untitled: CMakeFiles/untitled.dir/Iteration.cpp.o
untitled: CMakeFiles/untitled.dir/Population.cpp.o
untitled: CMakeFiles/untitled.dir/Chromosome.cpp.o
untitled: CMakeFiles/untitled.dir/python_test.cpp.o
untitled: CMakeFiles/untitled.dir/build.make
untitled: CMakeFiles/untitled.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable untitled"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/untitled.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/untitled.dir/build: untitled

.PHONY : CMakeFiles/untitled.dir/build

CMakeFiles/untitled.dir/requires: CMakeFiles/untitled.dir/main.cpp.o.requires
CMakeFiles/untitled.dir/requires: CMakeFiles/untitled.dir/Gene.cpp.o.requires
CMakeFiles/untitled.dir/requires: CMakeFiles/untitled.dir/Swarm.cpp.o.requires
CMakeFiles/untitled.dir/requires: CMakeFiles/untitled.dir/Iteration.cpp.o.requires
CMakeFiles/untitled.dir/requires: CMakeFiles/untitled.dir/Population.cpp.o.requires
CMakeFiles/untitled.dir/requires: CMakeFiles/untitled.dir/Chromosome.cpp.o.requires
CMakeFiles/untitled.dir/requires: CMakeFiles/untitled.dir/python_test.cpp.o.requires

.PHONY : CMakeFiles/untitled.dir/requires

CMakeFiles/untitled.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/untitled.dir/cmake_clean.cmake
.PHONY : CMakeFiles/untitled.dir/clean

CMakeFiles/untitled.dir/depend:
	cd /home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/biot/projects/szakdolgozat/Evolutionary_algorithm /home/biot/projects/szakdolgozat/Evolutionary_algorithm /home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug /home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug /home/biot/projects/szakdolgozat/Evolutionary_algorithm/cmake-build-debug/CMakeFiles/untitled.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/untitled.dir/depend

