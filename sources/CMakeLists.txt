cmake_minimum_required(VERSION 3.4)
project(Contig CXX)

#if (APPLE)
 #   include_directories(/opt/homebrew/Cellar/libomp/16.0.4/include)
  #  find_library(OPENMP_LIB libomp.dylib /opt/homebrew/Cellar/libomp/16.0.4/lib/)
   # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
#else()
 #   find_package(OpenMP REQUIRED)
  #  set(OPENMP_LIB OpenMP::OpenMP_CXX)
   # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif ()

#set(CMAKE_CXX_COMPILER g++-11)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wunused-parameter -O3")

# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")
# add the Sharg Parser to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../sharg-parser/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
find_package (sharg 1.0 REQUIRED)

######################### SOURCE ##############################

set(ULIB Utils)
set(FALIB FindAll)
set(CCLIB CodonCount)
set(CDLIB ContigDiff)

set(USRC ${ULIB}/directory.cpp ${ULIB}/fasta.cpp ${ULIB}/protein.cpp ${ULIB}/parser.cpp ${ULIB}/file.cpp)
set(FASRC ${FALIB}/find_all.cpp ${FALIB}/protein.cpp ${FALIB}/nucleic.cpp)
set(CCSRC ${CCLIB}/condo_count.cpp ${CCLIB}/parser.cpp)
set(CDSRC ${CDLIB}/contig_diff.cpp ${CDLIB}/parser.cpp ${CDLIB}/cent_percent.cpp ${CDLIB}/other_percent.cpp)

add_executable(Contig main.cpp ${USRC} ${FASRC} ${CCSRC} ${CDSRC})

#target_link_libraries(Contig LINK_PUBLIC ${OPENMP_LIB})
target_link_libraries(Contig seqan3::seqan3 sharg::sharg)

#EOF

