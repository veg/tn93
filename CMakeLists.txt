cmake_minimum_required(VERSION 3.0)
project(TN93)

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

set(CMAKE_CONFIGURATION_TYPES Release)

# SET VERSION FROM FILE
file(STRINGS "VERSION.txt" VERSION_NUMBER)
add_definitions(-DVERSION_NUMBER=\"${VERSION_NUMBER}\")
add_compile_options(-O3 -std=c++14 -funroll-loops)
add_link_options(-O3 -std=c++14 -funroll-loops)

include_directories(src/)

# List of sources for tn93
set(TN93_SOURCES src/TN93.cpp src/stringBuffer.cc src/tn93_shared.cc src/argparse.cpp)

# Add dynamically linked tn93 target
add_executable(tn93 ${TN93_SOURCES})
target_link_libraries(tn93 PRIVATE ${DEFAULT_LIBRARIES})

# Add statically linked tn93 target
add_executable(tn93_static.singlethread ${TN93_SOURCES})
target_link_options(tn93_static.singlethread PRIVATE -static-libgcc -static-libstdc++ -pthread -lm -lc)

# Add mixed statically linked tn93 target with static OpenMP
add_executable(tn93_mixed_openmp ${TN93_SOURCES})
target_link_options(tn93_mixed_openmp PRIVATE -static-libgcc -static-libstdc++ -pthread -lm -lc)

find_package(OpenMP REQUIRED)

if(OpenMP_FOUND)
    # Link OpenMP dynamically for dynamic target
    target_link_libraries(tn93 PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(tn93_mixed_openmp PRIVATE OpenMP::OpenMP_CXX)
endif()

# List of other executables and their sources
set(TN93_CLUSTER_SOURCES src/cluster.cpp src/stringBuffer.cc src/tn93_shared.cc src/argparse_cluster.cpp)
set(SEQCOVERAGE_SOURCES src/charfreqs.cpp src/stringBuffer.cc src/tn93_shared.cc src/argparse_cf.cpp)
set(READREDUCE_SOURCES src/read_reducer.cpp src/stringBuffer.cc src/tn93_shared.cc src/argparse_merge.cpp)
set(NUCFREQSFASTA_SOURCES src/nuc_freqs_from_fasta.cpp src/stringBuffer.cc src/tn93_shared.cc)
set(SHORTESTPATHTN93_SOURCES src/ShortestPathTN93.cpp src/stringBuffer.cc src/tn93_shared.cc)
set(SELECTREADS_SOURCES src/trim_reads.cpp src/stringBuffer.cc src/tn93_shared.cc src/argparse_trim.cpp)
set(VALIDATE_FASTA_SOURCES src/validate_fasta.cpp src/stringBuffer.cc src/tn93_shared.cc)
set(FASTA_DIGEST_SOURCES src/fasta_digest.cpp src/stringBuffer.cc src/tn93_shared.cc src/argparse_fasta_digest.cpp)
set(FASTA_DIFF_SOURCES src/fasta_diff.cpp src/stringBuffer.cc src/tn93_shared.cc src/argparse_fasta_diff.cpp)
set(SEQDIFF_SOURCES src/seqdiff.cpp src/stringBuffer.cc src/tn93_shared.cc src/argparse_seqdiff.cpp)

# Function to add dynamic executable
function(add_dynamic_executable target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE ${DEFAULT_LIBRARIES})
endfunction()

# Add other executables dynamically
add_dynamic_executable(tn93-cluster ${TN93_CLUSTER_SOURCES})
add_dynamic_executable(seqcoverage ${SEQCOVERAGE_SOURCES})
add_dynamic_executable(readreduce ${READREDUCE_SOURCES})
add_dynamic_executable(nucfreqsfasta ${NUCFREQSFASTA_SOURCES})
add_dynamic_executable(ShortestPathTN93 ${SHORTESTPATHTN93_SOURCES})
add_dynamic_executable(selectreads ${SELECTREADS_SOURCES})
add_dynamic_executable(validate_fasta ${VALIDATE_FASTA_SOURCES})
add_dynamic_executable(fasta_digest ${FASTA_DIGEST_SOURCES})
add_dynamic_executable(fasta_diff ${FASTA_DIFF_SOURCES})
add_dynamic_executable(seqdiff ${SEQDIFF_SOURCES})

# Install targets
install(TARGETS tn93 tn93_static.singlethread tn93_mixed_openmp tn93-cluster selectreads seqcoverage readreduce fasta_diff seqdiff fasta_digest ShortestPathTN93 nucfreqsfasta
    RUNTIME DESTINATION bin
    OPTIONAL
)

add_custom_target(TN93SP
    DEPENDS ShortestPathTN93
)

add_custom_target(NF
    DEPENDS nucfreqsfasta
)
