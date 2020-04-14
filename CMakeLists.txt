project(covid19-epi LANGUAGES C CUDA)
cmake_minimum_required(VERSION 3.13)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

set(COMMON_SOURCES COV_rand.c COV_rand.h common.h distance.c distance.h locale.c locale.h prop_distribution.c prop_distribution.h)

add_executable(covid19 covid19.c ${COMMON_SOURCES})
add_executable(covid19-gpu covid19.c ${COMMON_SOURCES} COV_gpu.cu)

target_compile_options(covid19 PRIVATE "-march=native" "-mtune=native")

target_compile_definitions(covid19 PRIVATE COV_GPU=0 USE_LOCALE_DISTANCE RMT_ENABLED=0)
target_compile_definitions(covid19-gpu PRIVATE COV_GPU=1 USE_LOCALE_DISTANCE RMT_ENABLED=0)

target_link_libraries(covid19 PRIVATE OpenMP::OpenMP_C Threads::Threads)
target_link_libraries(covid19-gpu PRIVATE OpenMP::OpenMP_C Threads::Threads)

if(UNIX)
    target_link_libraries(covid19 PRIVATE m dl)
    target_link_libraries(covid19-gpu PRIVATE m dl)
endif()

add_library(Remotery Remotery/lib/Remotery.c Remotery/lib/Remotery.h)
target_include_directories(Remotery INTERFACE Remotery/lib)
target_link_libraries(Remotery PRIVATE Threads::Threads m)

add_library(Remotery-gpu Remotery/lib/Remotery.c Remotery/lib/Remotery.h)
target_compile_definitions(Remotery-gpu PUBLIC RMT_USE_CUDA)
target_include_directories(Remotery-gpu INTERFACE Remotery/lib)
target_link_libraries(Remotery-gpu PRIVATE Threads::Threads m)

target_link_libraries(covid19 PRIVATE Remotery)
target_link_libraries(covid19-gpu PRIVATE Remotery-gpu)