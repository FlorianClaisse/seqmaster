cmake_minimum_required(VERSION 3.24)
project(Contig)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wunused-parameter")

######################### SOURCE ##############################

set(FLIB Foundation)
set(FALIB FindAll)

set(FSRC ${FLIB}/Directory.cpp ${FLIB}/FastaDecoder.cpp ${FLIB}/program_option.cpp FindAll/find_all.cpp FindAll/find_all.h)
set(FASRC ${FALIB}/find_all.cpp)

######################### LIBRARY ##############################

add_library(foundation ${FSRC})
add_library(findall ${FASRC})

######################### EXECUTABLE ##############################

add_executable(contig main.cpp)

######################### LINK LIBRARY ##############################

target_link_libraries(contig foundation)
target_link_libraries(contig findall)

#EOF

