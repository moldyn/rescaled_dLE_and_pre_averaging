#To "produce" a working binary, i.e., program, run the following two commands in the shell.
#cmake -H. -Bbuild
#cmake --build build -- -j3

cmake_minimum_required (VERSION 2.8)
project (LE_varygammavarym) 

# do we want static libraries?
# When STATIC_LINKING is set, cmake looks for libraries called *.a.
if(STATIC_LINKING)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif(STATIC_LINKING)



# Gnu compiler setup
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  message("using Gnu compiler")
  # use C++ 2011 standard
  set (CMAKE_CXX_FLAGS "-std=c++11")
  # optimization
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
  # parallelization
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  # warnings
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")                                                                          #-Wall does not give all error massages, so one needs the other -W... options
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmissing-declarations")
else()
  message(SEND_ERROR
    "currently no other compiler than the GCC is supported. sorry.")
endif()


#find_package (Boost 1.49 COMPONENTS program_options system filesystem REQUIRED)                                            #Flo uses boost for the program options
#include_directories (${Boost_INCLUDE_DIR})

find_package(Eigen3 REQUIRED)                                                                                               #Needed in script
include_directories (${EIGEN3_INCLUDE_DIR})


#add_subdirectory(subdirectory)


set(SRCS LE_varygammavarym.cpp)                                                                                           # Add here names of cource files

set(LIBS ${Boost_LIBRARIES})

## CUDA setup
#find_package(CUDA REQUIRED)                                                                                                #Stuff for CUDA -> graphic cards
#set(CUDA_PROPAGATE_HOST_FLAGS OFF)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}
#     "-arch=sm_30;-std=c++11;-O3;-DVERBOSE;-Xcompiler -fopenmp")
#set (LIBS ${LIBS} ${CUDA_LIBRARIES})
#include_directories (${CUDA_INCLUDE_DIRS})



set (PROGNAME LE_varygammavarym)


#cuda_add_executable(${PROGNAME} ${SRCS})                                                                                      # Signalword "cuda" nedded if CUDA code inclueded
add_executable(${PROGNAME} ${SRCS})
target_link_libraries(${PROGNAME} ${LIBS})


# set -static, when STATIC_LINKING is TRUE and set LINK_SEARCH_END_STATIC
# to remove the additional -bdynamic from the linker line.
if(STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
    set_target_properties(${PROGNAME} PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(STATIC_LINKING)


install(TARGETS ${PROGNAME} RUNTIME DESTINATION .)

# create source-doc with doxygen
add_custom_target(doc
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  COMMAND doxygen doxygen.config)

