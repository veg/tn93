# Docker image for a tn93 development environment
FROM oraclelinux:8

# Set up environment and install dependencies
RUN yum -y update && \
    yum install -y cmake gcc-c++ gcc-toolset-10 git make oracle-epel-release-el8 && \
    rm -f CMakeCache.txt && \
    scl enable gcc-toolset-10 bash

# To compile tn93 within the development environment:
#   cmake .
#   make
