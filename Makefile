# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /apps/Hebbe7/software/Compiler/GCCcore/8.3.0/CMake/3.15.3/bin/cmake

# The command to remove a file.
RM = /apps/Hebbe7/software/Compiler/GCCcore/8.3.0/CMake/3.15.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cephyr/users/ssoheil/Hebbe/work_dir/covid19-repo-v5-branched

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cephyr/users/ssoheil/Hebbe/work_dir/covid19-repo-v5-branched

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/apps/Hebbe7/software/Compiler/GCCcore/8.3.0/CMake/3.15.3/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/apps/Hebbe7/software/Compiler/GCCcore/8.3.0/CMake/3.15.3/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /cephyr/users/ssoheil/Hebbe/work_dir/covid19-repo-v5-branched/CMakeFiles /cephyr/users/ssoheil/Hebbe/work_dir/covid19-repo-v5-branched/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /cephyr/users/ssoheil/Hebbe/work_dir/covid19-repo-v5-branched/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named covid19

# Build rule for target.
covid19: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 covid19
.PHONY : covid19

# fast build rule for target.
covid19/fast:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/build
.PHONY : covid19/fast

COV_rand.o: COV_rand.c.o

.PHONY : COV_rand.o

# target to build an object file
COV_rand.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/COV_rand.c.o
.PHONY : COV_rand.c.o

COV_rand.i: COV_rand.c.i

.PHONY : COV_rand.i

# target to preprocess a source file
COV_rand.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/COV_rand.c.i
.PHONY : COV_rand.c.i

COV_rand.s: COV_rand.c.s

.PHONY : COV_rand.s

# target to generate assembly for a file
COV_rand.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/COV_rand.c.s
.PHONY : COV_rand.c.s

age_dist.o: age_dist.c.o

.PHONY : age_dist.o

# target to build an object file
age_dist.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/age_dist.c.o
.PHONY : age_dist.c.o

age_dist.i: age_dist.c.i

.PHONY : age_dist.i

# target to preprocess a source file
age_dist.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/age_dist.c.i
.PHONY : age_dist.c.i

age_dist.s: age_dist.c.s

.PHONY : age_dist.s

# target to generate assembly for a file
age_dist.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/age_dist.c.s
.PHONY : age_dist.c.s

calc_infect.o: calc_infect.c.o

.PHONY : calc_infect.o

# target to build an object file
calc_infect.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/calc_infect.c.o
.PHONY : calc_infect.c.o

calc_infect.i: calc_infect.c.i

.PHONY : calc_infect.i

# target to preprocess a source file
calc_infect.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/calc_infect.c.i
.PHONY : calc_infect.c.i

calc_infect.s: calc_infect.c.s

.PHONY : calc_infect.s

# target to generate assembly for a file
calc_infect.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/calc_infect.c.s
.PHONY : calc_infect.c.s

calc_kappa.o: calc_kappa.c.o

.PHONY : calc_kappa.o

# target to build an object file
calc_kappa.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/calc_kappa.c.o
.PHONY : calc_kappa.c.o

calc_kappa.i: calc_kappa.c.i

.PHONY : calc_kappa.i

# target to preprocess a source file
calc_kappa.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/calc_kappa.c.i
.PHONY : calc_kappa.c.i

calc_kappa.s: calc_kappa.c.s

.PHONY : calc_kappa.s

# target to generate assembly for a file
calc_kappa.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/calc_kappa.c.s
.PHONY : calc_kappa.c.s

city_lat_long.o: city_lat_long.c.o

.PHONY : city_lat_long.o

# target to build an object file
city_lat_long.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/city_lat_long.c.o
.PHONY : city_lat_long.c.o

city_lat_long.i: city_lat_long.c.i

.PHONY : city_lat_long.i

# target to preprocess a source file
city_lat_long.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/city_lat_long.c.i
.PHONY : city_lat_long.c.i

city_lat_long.s: city_lat_long.c.s

.PHONY : city_lat_long.s

# target to generate assembly for a file
city_lat_long.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/city_lat_long.c.s
.PHONY : city_lat_long.c.s

covid19.o: covid19.c.o

.PHONY : covid19.o

# target to build an object file
covid19.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/covid19.c.o
.PHONY : covid19.c.o

covid19.i: covid19.c.i

.PHONY : covid19.i

# target to preprocess a source file
covid19.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/covid19.c.i
.PHONY : covid19.c.i

covid19.s: covid19.c.s

.PHONY : covid19.s

# target to generate assembly for a file
covid19.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/covid19.c.s
.PHONY : covid19.c.s

death.o: death.c.o

.PHONY : death.o

# target to build an object file
death.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/death.c.o
.PHONY : death.c.o

death.i: death.c.i

.PHONY : death.i

# target to preprocess a source file
death.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/death.c.i
.PHONY : death.c.i

death.s: death.c.s

.PHONY : death.s

# target to generate assembly for a file
death.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/death.c.s
.PHONY : death.c.s

distance.o: distance.c.o

.PHONY : distance.o

# target to build an object file
distance.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/distance.c.o
.PHONY : distance.c.o

distance.i: distance.c.i

.PHONY : distance.i

# target to preprocess a source file
distance.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/distance.c.i
.PHONY : distance.c.i

distance.s: distance.c.s

.PHONY : distance.s

# target to generate assembly for a file
distance.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/distance.c.s
.PHONY : distance.c.s

hospital.o: hospital.c.o

.PHONY : hospital.o

# target to build an object file
hospital.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/hospital.c.o
.PHONY : hospital.c.o

hospital.i: hospital.c.i

.PHONY : hospital.i

# target to preprocess a source file
hospital.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/hospital.c.i
.PHONY : hospital.c.i

hospital.s: hospital.c.s

.PHONY : hospital.s

# target to generate assembly for a file
hospital.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/hospital.c.s
.PHONY : hospital.c.s

household_lat_long.o: household_lat_long.c.o

.PHONY : household_lat_long.o

# target to build an object file
household_lat_long.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/household_lat_long.c.o
.PHONY : household_lat_long.c.o

household_lat_long.i: household_lat_long.c.i

.PHONY : household_lat_long.i

# target to preprocess a source file
household_lat_long.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/household_lat_long.c.i
.PHONY : household_lat_long.c.i

household_lat_long.s: household_lat_long.c.s

.PHONY : household_lat_long.s

# target to generate assembly for a file
household_lat_long.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/household_lat_long.c.s
.PHONY : household_lat_long.c.s

initialize_infections.o: initialize_infections.c.o

.PHONY : initialize_infections.o

# target to build an object file
initialize_infections.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/initialize_infections.c.o
.PHONY : initialize_infections.c.o

initialize_infections.i: initialize_infections.c.i

.PHONY : initialize_infections.i

# target to preprocess a source file
initialize_infections.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/initialize_infections.c.i
.PHONY : initialize_infections.c.i

initialize_infections.s: initialize_infections.c.s

.PHONY : initialize_infections.s

# target to generate assembly for a file
initialize_infections.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/initialize_infections.c.s
.PHONY : initialize_infections.c.s

job_dist.o: job_dist.c.o

.PHONY : job_dist.o

# target to build an object file
job_dist.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/job_dist.c.o
.PHONY : job_dist.c.o

job_dist.i: job_dist.c.i

.PHONY : job_dist.i

# target to preprocess a source file
job_dist.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/job_dist.c.i
.PHONY : job_dist.c.i

job_dist.s: job_dist.c.s

.PHONY : job_dist.s

# target to generate assembly for a file
job_dist.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/job_dist.c.s
.PHONY : job_dist.c.s

locale.o: locale.c.o

.PHONY : locale.o

# target to build an object file
locale.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/locale.c.o
.PHONY : locale.c.o

locale.i: locale.c.i

.PHONY : locale.i

# target to preprocess a source file
locale.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/locale.c.i
.PHONY : locale.c.i

locale.s: locale.c.s

.PHONY : locale.s

# target to generate assembly for a file
locale.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/locale.c.s
.PHONY : locale.c.s

prop_distribution.o: prop_distribution.c.o

.PHONY : prop_distribution.o

# target to build an object file
prop_distribution.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/prop_distribution.c.o
.PHONY : prop_distribution.c.o

prop_distribution.i: prop_distribution.c.i

.PHONY : prop_distribution.i

# target to preprocess a source file
prop_distribution.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/prop_distribution.c.i
.PHONY : prop_distribution.c.i

prop_distribution.s: prop_distribution.c.s

.PHONY : prop_distribution.s

# target to generate assembly for a file
prop_distribution.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/prop_distribution.c.s
.PHONY : prop_distribution.c.s

segment_population.o: segment_population.c.o

.PHONY : segment_population.o

# target to build an object file
segment_population.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/segment_population.c.o
.PHONY : segment_population.c.o

segment_population.i: segment_population.c.i

.PHONY : segment_population.i

# target to preprocess a source file
segment_population.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/segment_population.c.i
.PHONY : segment_population.c.i

segment_population.s: segment_population.c.s

.PHONY : segment_population.s

# target to generate assembly for a file
segment_population.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/segment_population.c.s
.PHONY : segment_population.c.s

workplace_dist.o: workplace_dist.c.o

.PHONY : workplace_dist.o

# target to build an object file
workplace_dist.c.o:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/workplace_dist.c.o
.PHONY : workplace_dist.c.o

workplace_dist.i: workplace_dist.c.i

.PHONY : workplace_dist.i

# target to preprocess a source file
workplace_dist.c.i:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/workplace_dist.c.i
.PHONY : workplace_dist.c.i

workplace_dist.s: workplace_dist.c.s

.PHONY : workplace_dist.s

# target to generate assembly for a file
workplace_dist.c.s:
	$(MAKE) -f CMakeFiles/covid19.dir/build.make CMakeFiles/covid19.dir/workplace_dist.c.s
.PHONY : workplace_dist.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... covid19"
	@echo "... edit_cache"
	@echo "... COV_rand.o"
	@echo "... COV_rand.i"
	@echo "... COV_rand.s"
	@echo "... age_dist.o"
	@echo "... age_dist.i"
	@echo "... age_dist.s"
	@echo "... calc_infect.o"
	@echo "... calc_infect.i"
	@echo "... calc_infect.s"
	@echo "... calc_kappa.o"
	@echo "... calc_kappa.i"
	@echo "... calc_kappa.s"
	@echo "... city_lat_long.o"
	@echo "... city_lat_long.i"
	@echo "... city_lat_long.s"
	@echo "... covid19.o"
	@echo "... covid19.i"
	@echo "... covid19.s"
	@echo "... death.o"
	@echo "... death.i"
	@echo "... death.s"
	@echo "... distance.o"
	@echo "... distance.i"
	@echo "... distance.s"
	@echo "... hospital.o"
	@echo "... hospital.i"
	@echo "... hospital.s"
	@echo "... household_lat_long.o"
	@echo "... household_lat_long.i"
	@echo "... household_lat_long.s"
	@echo "... initialize_infections.o"
	@echo "... initialize_infections.i"
	@echo "... initialize_infections.s"
	@echo "... job_dist.o"
	@echo "... job_dist.i"
	@echo "... job_dist.s"
	@echo "... locale.o"
	@echo "... locale.i"
	@echo "... locale.s"
	@echo "... prop_distribution.o"
	@echo "... prop_distribution.i"
	@echo "... prop_distribution.s"
	@echo "... segment_population.o"
	@echo "... segment_population.i"
	@echo "... segment_population.s"
	@echo "... workplace_dist.o"
	@echo "... workplace_dist.i"
	@echo "... workplace_dist.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

