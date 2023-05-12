cmake_minimum_required(VERSION 3.10)
project(Contig)

if (APPLE)
    include_directories(/opt/homebrew/Cellar/libomp/16.0.3/include)
    find_library(OPENMP_LIB libomp.dylib /opt/homebrew/Cellar/libomp/16.0.3/lib/)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
else()
    find_package(OpenMP REQUIRED)
    set(OPENMP_LIB OpenMP::OpenMP_CXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wunused-parameter -O3")

######################### SOURCE ##############################

set(FLIB Foundation)
set(FALIB FindAll)
set(CCLIB CodonCount)
set(CDLIB ContigDiff)

set(FSRC ${FLIB}/directory.cpp ${FLIB}/fasta.cpp ${FLIB}/program_option.cpp)
set(FASRC ${FALIB}/find_all.cpp)
set(CCSRC ${CCLIB}/condo_count.cpp)
set(CDSRC ${CDLIB}/contig_diff.cpp ${CDLIB}/search_algo.cpp)

add_executable(Contig main.cpp ${FSRC} ${FASRC} ${CCSRC} ${CDSRC})

target_link_libraries(Contig LINK_PUBLIC ${OPENMP_LIB})

#EOF

