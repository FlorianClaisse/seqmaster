cmake_minimum_required(VERSION 3.4)
project(Contig CXX)

# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")
# add the Sharg Parser to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../sharg-parser/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
find_package (sharg 1.0 REQUIRED)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/csv-parser)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/csv-parser/include)

######################### SOURCE ##############################

set(ULIB Utils)
set(FALIB FindAll)
set(CCLIB CodonCount)
set(CDLIB ContigDiff)
set(GMLIB GeneMut)

set(USRC ${ULIB}/directory.cpp ${ULIB}/protein.cpp ${ULIB}/parser.cpp ${ULIB}/file.cpp)
set(FASRC ${FALIB}/find_all.cpp)
set(CCSRC ${CCLIB}/condo_count.cpp)
set(CDSRC ${CDLIB}/contig_diff.cpp)
set(GMSRC ${GMLIB}/gene_mut.cpp)

add_executable(Contig main.cpp ${USRC} ${FASRC} ${CCSRC} ${CDSRC} ${GMSRC})

target_link_libraries(Contig seqan3::seqan3 sharg::sharg csv)

#EOF

