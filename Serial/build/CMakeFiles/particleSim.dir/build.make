# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

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
CMAKE_SOURCE_DIR = /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/build

# Include any dependencies generated for this target.
include CMakeFiles/particleSim.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/particleSim.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/particleSim.dir/flags.make

CMakeFiles/particleSim.dir/main.cpp.o: CMakeFiles/particleSim.dir/flags.make
CMakeFiles/particleSim.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/particleSim.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/particleSim.dir/main.cpp.o -c /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/main.cpp

CMakeFiles/particleSim.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/particleSim.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/main.cpp > CMakeFiles/particleSim.dir/main.cpp.i

CMakeFiles/particleSim.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/particleSim.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/main.cpp -o CMakeFiles/particleSim.dir/main.cpp.s

CMakeFiles/particleSim.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/particleSim.dir/main.cpp.o.requires

CMakeFiles/particleSim.dir/main.cpp.o.provides: CMakeFiles/particleSim.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/particleSim.dir/build.make CMakeFiles/particleSim.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/particleSim.dir/main.cpp.o.provides

CMakeFiles/particleSim.dir/main.cpp.o.provides.build: CMakeFiles/particleSim.dir/main.cpp.o

CMakeFiles/particleSim.dir/particle.cpp.o: CMakeFiles/particleSim.dir/flags.make
CMakeFiles/particleSim.dir/particle.cpp.o: ../particle.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/particleSim.dir/particle.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/particleSim.dir/particle.cpp.o -c /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/particle.cpp

CMakeFiles/particleSim.dir/particle.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/particleSim.dir/particle.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/particle.cpp > CMakeFiles/particleSim.dir/particle.cpp.i

CMakeFiles/particleSim.dir/particle.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/particleSim.dir/particle.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/particle.cpp -o CMakeFiles/particleSim.dir/particle.cpp.s

CMakeFiles/particleSim.dir/particle.cpp.o.requires:
.PHONY : CMakeFiles/particleSim.dir/particle.cpp.o.requires

CMakeFiles/particleSim.dir/particle.cpp.o.provides: CMakeFiles/particleSim.dir/particle.cpp.o.requires
	$(MAKE) -f CMakeFiles/particleSim.dir/build.make CMakeFiles/particleSim.dir/particle.cpp.o.provides.build
.PHONY : CMakeFiles/particleSim.dir/particle.cpp.o.provides

CMakeFiles/particleSim.dir/particle.cpp.o.provides.build: CMakeFiles/particleSim.dir/particle.cpp.o

CMakeFiles/particleSim.dir/particles.cpp.o: CMakeFiles/particleSim.dir/flags.make
CMakeFiles/particleSim.dir/particles.cpp.o: ../particles.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/particleSim.dir/particles.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/particleSim.dir/particles.cpp.o -c /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/particles.cpp

CMakeFiles/particleSim.dir/particles.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/particleSim.dir/particles.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/particles.cpp > CMakeFiles/particleSim.dir/particles.cpp.i

CMakeFiles/particleSim.dir/particles.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/particleSim.dir/particles.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/particles.cpp -o CMakeFiles/particleSim.dir/particles.cpp.s

CMakeFiles/particleSim.dir/particles.cpp.o.requires:
.PHONY : CMakeFiles/particleSim.dir/particles.cpp.o.requires

CMakeFiles/particleSim.dir/particles.cpp.o.provides: CMakeFiles/particleSim.dir/particles.cpp.o.requires
	$(MAKE) -f CMakeFiles/particleSim.dir/build.make CMakeFiles/particleSim.dir/particles.cpp.o.provides.build
.PHONY : CMakeFiles/particleSim.dir/particles.cpp.o.provides

CMakeFiles/particleSim.dir/particles.cpp.o.provides.build: CMakeFiles/particleSim.dir/particles.cpp.o

# Object files for target particleSim
particleSim_OBJECTS = \
"CMakeFiles/particleSim.dir/main.cpp.o" \
"CMakeFiles/particleSim.dir/particle.cpp.o" \
"CMakeFiles/particleSim.dir/particles.cpp.o"

# External object files for target particleSim
particleSim_EXTERNAL_OBJECTS =

particleSim: CMakeFiles/particleSim.dir/main.cpp.o
particleSim: CMakeFiles/particleSim.dir/particle.cpp.o
particleSim: CMakeFiles/particleSim.dir/particles.cpp.o
particleSim: CMakeFiles/particleSim.dir/build.make
particleSim: CMakeFiles/particleSim.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable particleSim"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/particleSim.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/particleSim.dir/build: particleSim
.PHONY : CMakeFiles/particleSim.dir/build

CMakeFiles/particleSim.dir/requires: CMakeFiles/particleSim.dir/main.cpp.o.requires
CMakeFiles/particleSim.dir/requires: CMakeFiles/particleSim.dir/particle.cpp.o.requires
CMakeFiles/particleSim.dir/requires: CMakeFiles/particleSim.dir/particles.cpp.o.requires
.PHONY : CMakeFiles/particleSim.dir/requires

CMakeFiles/particleSim.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/particleSim.dir/cmake_clean.cmake
.PHONY : CMakeFiles/particleSim.dir/clean

CMakeFiles/particleSim.dir/depend:
	cd /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/build /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/build /home/Student/s4638706/COSC3500/COSC3500ParticleSim/Serial/build/CMakeFiles/particleSim.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/particleSim.dir/depend

