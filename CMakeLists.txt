cmake_minimum_required(VERSION 3.10)
project(Contig)

if (APPLE)
    include_directories(/opt/homebrew/Cellar/libomp/16.0.3/include)
    include_directories(/usr/local/Cellar/libomp/16.0.3/include)
    find_library(OPENMP_LIB libomp.dylib /opt/homebrew/Cellar/libomp/16.0.3/lib/)
    find_library(OPENMP_LIB libomp.dylib /usr/local/Cellar/libomp/16.0.3/lib)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
else()
    find_package(OpenMP REQUIRED)
    set(OPENMP_LIB OpenMP::OpenMP_CXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wunused-parameter -O3")

######################### SOURCE ##############################

set(ULIB Utils)
set(FALIB FindAll)
set(CCLIB CodonCount)
set(CDLIB ContigDiff)

set(USRC ${ULIB}/directory.cpp ${ULIB}/fasta.cpp)
set(FASRC ${FALIB}/find_all.cpp ${FALIB}/parser.cpp)
set(CCSRC ${CCLIB}/condo_count.cpp ${CCLIB}/parser.cpp)
set(CDSRC ${CDLIB}/contig_diff.cpp ${CDLIB}/parser.cpp)

add_executable(Contig main.cpp ${USRC} ${FASRC} ${CCSRC} ${CDSRC})

target_link_libraries(Contig LINK_PUBLIC ${OPENMP_LIB})

#EOF

