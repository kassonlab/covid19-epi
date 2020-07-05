project(covid19-epi LANGUAGES C CXX)
cmake_minimum_required(VERSION 3.10)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -xHost -ipo" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -xHost -ipo" CACHE STRING "" FORCE)
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -xHost -ipo" CACHE STRING "" FORCE)
    set(CMAKE_CUDA_FLAGS_RELEASE_INIT "${CMAKE_CUDA_FLAGS_RELEASE_INIT} --compiler-options=-xHost --compiler-options=-O3" CACHE STRING "" FORCE)
    set(MATH_LIBS "imf")
    set(ENV{CUDAHOSTCXX} "icpc")
else()
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -mtune=native" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native" CACHE STRING "" FORCE)
    set(CMAKE_CUDA_FLAGS_RELEASE_INIT "${CMAKE_CUDA_FLAGS_RELEASE_INIT} --compiler-options=-march=native --compiler-options=-mtune=native --compiler-options=-O3" CACHE STRING "" FORCE)
    set(MATH_LIBS "m")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo." FORCE)
endif()

include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    if(NOT GPU_ARCHITECTURE)
        set(GPU_ARCHITECTURE sm_70)
    endif()
    message("GPU architecture selected: ${GPU_ARCHITECTURE}")
endif()

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

set(COMMON_SOURCES COV_rand.c COV_rand.h common.h distance.c distance.h interventions.c infections.c infections.h locale.c locale.h prop_distribution.c prop_distribution.h)

add_executable(covid19 covid19.c ${COMMON_SOURCES})
target_compile_options(covid19 PRIVATE ${ARCH_FLAGS} ${TUNE_FLAGS})
target_compile_definitions(covid19 PRIVATE COV_GPU=0 USE_LOCALE_DISTANCE)
target_link_libraries(covid19 PRIVATE OpenMP::OpenMP_C Threads::Threads)

if(UNIX)
    target_link_libraries(covid19 PRIVATE ${MATH_LIBS} dl)
endif()

if(CMAKE_CUDA_COMPILER)
    add_executable(covid19-gpu covid19.c ${COMMON_SOURCES} COV_gpu.cu)

    target_compile_options(covid19-gpu PRIVATE
        "$<$<COMPILE_LANGUAGE:CUDA>:--gpu-architecture=${GPU_ARCHITECTURE}>")

    target_compile_definitions(covid19-gpu PRIVATE COV_GPU=1 USE_LOCALE_DISTANCE)
    target_link_libraries(covid19-gpu PRIVATE OpenMP::OpenMP_C Threads::Threads)

    if(UNIX)
        target_link_libraries(covid19-gpu PRIVATE ${MATH_LIBS} dl)
    endif()
endif()
