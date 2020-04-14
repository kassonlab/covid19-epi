project(covid19-epi LANGUAGES C CXX)
cmake_minimum_required(VERSION 3.13)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

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

set(COMMON_SOURCES COV_rand.c COV_rand.h common.h distance.c distance.h locale.c locale.h prop_distribution.c prop_distribution.h)

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    set(ARCH_FLAGS "-xHost")
    set(TUNE_FLAGS "")
    set(MATH_LIBS "imf")
else()
    set(ARCH_FLAGS "-march=native")
    set(TUNE_FLAGS "-mtune=native")
    set(MATH_LIBS "m")
endif()

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
        "$<$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>:${ARCH_FLAGS}>"
        "$<$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>:${TUNE_FLAGS}>"
        "$<$<COMPILE_LANGUAGE:CUDA>:-ccbin=${CMAKE_CXX_COMPILER}>"
        "$<$<COMPILE_LANGUAGE:CUDA>:--gpu-architecture=${GPU_ARCHITECTURE}>"
        "$<$<COMPILE_LANGUAGE:CUDA>:--compiler-options=${ARCH_FLAGS}>"
        "$<$<COMPILE_LANGUAGE:CUDA>:--compiler-options=${TUNE_FLAGS}>")

    target_compile_definitions(covid19-gpu PRIVATE COV_GPU=1 USE_LOCALE_DISTANCE)
    target_link_libraries(covid19-gpu PRIVATE OpenMP::OpenMP_C Threads::Threads)

    if(UNIX)
        target_link_libraries(covid19-gpu PRIVATE ${MATH_LIBS} dl)
    endif()
endif()
