cmake_minimum_required(VERSION 3.0)
project(TN93)

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

set(CMAKE_CONFIGURATION_TYPES Release)

# SET VERSION FROM FILE
file (STRINGS "VERSION.txt" VERSION_NUMBER)
add_definitions (-DVERSION_NUMBER=\"${VERSION_NUMBER}\")
add_compile_options (-O3 -std=c++14 -march=native -funroll-loops )
add_link_options(-O3 -std=c++14 -march=native -funroll-loops)

include_directories(
    src/
)

add_executable(
    tn93
    src/TN93.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse.cpp
)

add_executable(
    tn93-cluster
    src/cluster.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse_cluster.cpp
)

add_executable(
    seqcoverage
    src/charfreqs.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse_cf.cpp
)

add_executable(
    readreduce
    src/read_reducer.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse_merge.cpp
)

add_executable(
    nucfreqsfasta
    src/nuc_freqs_from_fasta.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
)

add_executable(
    ShortestPathTN93
    src/ShortestPathTN93.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
)

add_executable(
    selectreads
    src/trim_reads.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse_trim.cpp
)

add_executable(
    validate_fasta
    src/validate_fasta.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
)

add_executable(
    fasta_diff
    src/fasta_diff.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse_fasta_diff.cpp
)

add_executable(
    seqdiff
    src/seqdiff.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse_seqdiff.cpp
)


#set(DEFAULT_COMPILE_OPTIONS -O3 -std=c++14 -funroll-loops )
#set(DEFAULT_COMPILE_OPTIONS -std=c++14 -g -fsanitize=address -O0)
#set(DEFAULT_LINK_OPTIONS ${DEFAULT_COMPILE_OPTIONS})



#add_definitions(${DEFAULT_COMPILE_OPTIONS})

find_package(OpenMP)

if(${OPENMP_FOUND})
   add_definitions(${OpenMP_CXX_FLAGS})
   target_link_libraries(tn93 PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(tn93-cluster PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(selectreads PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(seqcoverage PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(readreduce PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(fasta_diff PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(ShortestPathTN93 PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(seqdiff PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
endif(${OPENMP_FOUND})



install(
    TARGETS tn93
	RUNTIME DESTINATION bin
	OPTIONAL
)

install(
    TARGETS tn93-cluster
	RUNTIME DESTINATION bin
	OPTIONAL
)

install(
    TARGETS selectreads
	RUNTIME DESTINATION bin
	OPTIONAL
)

install(
    TARGETS seqcoverage
	RUNTIME DESTINATION bin
	OPTIONAL
)

install(
    TARGETS readreduce
	RUNTIME DESTINATION bin
	OPTIONAL
)

install(
    TARGETS fasta_diff
	RUNTIME DESTINATION bin
	OPTIONAL
)

install(
    TARGETS seqdiff
	RUNTIME DESTINATION bin
	OPTIONAL
)

#add_custom_target(TN93
#  DEPENDS tn93
#)


install(
    TARGETS ShortestPathTN93
	RUNTIME DESTINATION bin
	OPTIONAL
)

add_custom_target(TN93SP
	DEPENDS ShortestPathTN93
)


install(
    TARGETS nucfreqsfasta
	RUNTIME DESTINATION bin
	OPTIONAL
)

add_custom_target(NF
	DEPENDS nucfreqsfasta
)
