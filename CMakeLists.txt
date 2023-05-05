cmake_minimum_required(VERSION 3.10)
project(Contig)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wunused-parameter")

######################### SOURCE ##############################

set(FLIB Foundation)
set(FALIB FindAll)
set(CCLIB CodonCount)

set(FSRC ${FLIB}/directory.cpp ${FLIB}/fasta.cpp ${FLIB}/program_option.cpp)
set(FASRC ${FALIB}/find_all.cpp)
set(CCSRC ${CCLIB}/condo_count.cpp)

add_executable(Contig main.cpp ${FSRC} ${FASRC} ${CCSRC})

#EOF

