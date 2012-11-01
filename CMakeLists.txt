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
    TN93dist
    src/TN93.cpp
    src/stringBuffer.cc
)

add_executable(
    ShortestPathTN93
    src/ShortestPathTN93.cpp
    src/stringBuffer.cc
)

set(DEFAULT_COMPILE_FLAGS "-O3")
set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")

find_package(OpenMP)

if (${OPENMP_FOUND})
	set_target_properties(
		TN93dist
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
endif (${OPENMP_FOUND})

install(
    TARGETS TN93dist
	RUNTIME DESTINATION bin
	OPTIONAL
)

add_custom_target(TN93 
	DEPENDS TN93dist
)


install(
    TARGETS ShortestPathTN93
	RUNTIME DESTINATION bin
	OPTIONAL
)

add_custom_target(TN93SP
	DEPENDS ShortestPathTN93
)

