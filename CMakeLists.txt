cmake_minimum_required(VERSION 3.24)
project(Contig)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wunused-parameter")

######################### SOURCE ##############################

set(FLIB Foundation)
set(FALIB FindAll)

set(FSRC ${FLIB}/directory.cpp ${FLIB}/Headers/directory.h ${FLIB}/fasta.cpp ${FLIB}/Headers/fasta.h ${FLIB}/program_option.cpp ${FLIB}/Headers/program_option.h)
set(FASRC ${FALIB}/find_all.cpp ${FALIB}/find_all.h)

add_executable(Contig main.cpp ${FSRC} ${FASRC})

#EOF

