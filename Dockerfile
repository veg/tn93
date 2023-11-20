# Docker image for a tn93 development environment
FROM oraclelinux:8

# Set up environment and install dependencies
RUN yum -y update && \
    yum install -y cmake gcc-c++ gcc-toolset-10 git make oracle-epel-release-el8 && \
    echo 'source /opt/rh/gcc-toolset-10/enable' > ~/.bashrc && \
    source ~/.bashrc

# To compile tn93 within the development environment:
#   cmake .
#   make
