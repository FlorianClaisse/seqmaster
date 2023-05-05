cmake_minimum_required(VERSION 3.10)
project(Contig)

#if(NOT TARGET OpenMP::OpenMP_CXX)
 #   find_package(Threads REQUIRED)
  #  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    #set_property(TARGET OpenMP::OpenMP_CXX
     #       PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    #set_property(TARGET OpenMP::OpenMP_CXX
     #       PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

#endif()

include_directories(/opt/homebrew/Cellar/libomp/16.0.3/include)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
#find_library(OPENMP_LIB )
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -Xpreprocessor -fopenmp")
find_library(OPENMP_LIB libomp.dylib /opt/homebrew/Cellar/libomp/16.0.3/lib/)
######################### SOURCE ##############################

set(FLIB Foundation)
set(FALIB FindAll)
set(CCLIB CodonCount)

set(FSRC ${FLIB}/directory.cpp ${FLIB}/fasta.cpp ${FLIB}/program_option.cpp)
set(FASRC ${FALIB}/find_all.cpp)
set(CCSRC ${CCLIB}/condo_count.cpp)

add_executable(Contig main.cpp ${FSRC} ${FASRC} ${CCSRC})

target_link_libraries(Contig LINK_PUBLIC ${OPENMP_LIB})

#EOF

