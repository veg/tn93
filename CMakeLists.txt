cmake_minimum_required(VERSION 2.8)
project(BinMix)

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
    BinMix
    src/BinomialMixture.cpp
)


install(
    TARGETS BinMix
	RUNTIME DESTINATION bin
	OPTIONAL
)

