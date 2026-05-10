cmake_minimum_required(VERSION 3.5)
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
add_compile_options (-O3 -std=c++14  -funroll-loops )
add_link_options(-O3 -std=c++14 -funroll-loops)

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
    fasta_digest
    src/fasta_digest.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse_fasta_digest.cpp
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

install(
    TARGETS fasta_digest
	RUNTIME DESTINATION bin
	OPTIONAL
)

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


enable_testing()

set(TN93_BIN ${CMAKE_CURRENT_BINARY_DIR}/tn93)
set(REF_TN93_BIN /usr/local/bin/tn93)
set(TEST_RUNNER ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_runner.py)
set(DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(GROUND_TRUTH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests/ground_truth)

find_package(PythonInterp 3 REQUIRED)

# Test 1: Small dataset, default settings
add_test(NAME small_default COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/test.fas --name small_default --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 2: Small dataset, resolve ambiguities
add_test(NAME small_resolve COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/test.fas --args "-a resolve" --name small_resolve --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 3: Small dataset, average ambiguities
add_test(NAME small_average COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/test.fas --args "-a average" --name small_average --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 4: Small dataset, Hamming skip
add_test(NAME small_hamming COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/test.fas --args "-H -t 0.05" --name small_hamming --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 5: Small dataset, custom overlap
add_test(NAME small_overlap COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/test.fas --args "-l 50" --name small_overlap --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 6: Subset dataset, default
add_test(NAME subset_default COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/subset.fasta --name subset_default --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 7: Subset dataset, Hamming skip
add_test(NAME subset_hamming COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/subset.fasta --args "-H -t 0.015" --name subset_hamming --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 8: Small dataset, csvn format
add_test(NAME small_csvn COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/test.fas --args "-f csvn" --name small_csvn --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 9: Small dataset, gapmm
add_test(NAME small_gapmm COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/test.fas --args "-a gapmm" --name small_gapmm --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 10: Large dataset, default
add_test(NAME large_default COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/large.fasta --name large_default --ground-truth-dir ${GROUND_TRUTH_DIR})

# Test 11: Large dataset, Hamming skip
add_test(NAME large_hamming COMMAND ${PYTHON_EXECUTABLE} ${TEST_RUNNER} 
    --binary ${TN93_BIN} --ref-binary ${REF_TN93_BIN} --input ${DATA_DIR}/large.fasta --args "-H -t 0.015" --name large_hamming --ground-truth-dir ${GROUND_TRUTH_DIR})
