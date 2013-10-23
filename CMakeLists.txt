cmake_minimum_required(VERSION 2.8)
project(TN93)

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)


set(CMAKE_CONFIGURATION_TYPES Release)

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

set(DEFAULT_COMPILE_FLAGS "-O3")
set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")

find_package(OpenMP)

if (${OPENMP_FOUND})
	set_target_properties(
		tn93
		PROPERTIES
		COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
		LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
	)
	set_target_properties(
		selectreads
		PROPERTIES
		COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
		LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
	)
	set_target_properties(
		seqcoverage
		PROPERTIES
		COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
		LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
	)
	set_target_properties(
		ShortestPathTN93
		PROPERTIES
		COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
		LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
    )
	set_target_properties(
		readreduce
		PROPERTIES
		COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
		LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
    )
endif (${OPENMP_FOUND})

install(
    TARGETS tn93
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

add_custom_target(TN93 
	DEPENDS tn93
)


install(
    TARGETS ShortestPathTN93
	RUNTIME DESTINATION bin
	OPTIONAL
)

add_custom_target(TN93SP
	DEPENDS ShortestPathTN93
)


add_custom_target(SEQCOVERAGE
	DEPENDS seqcoverage
)

install(
    TARGETS nucfreqsfasta
	RUNTIME DESTINATION bin
	OPTIONAL
)

add_custom_target(NF
	DEPENDS nucfreqsfasta
)
