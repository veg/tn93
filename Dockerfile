# Docker image for a tn93 development environment
FROM oraclelinux:8

# Set up environment and install dependencies
RUN yum -y check-update && yum -y update && \
    yum install -y cmake gcc-c++ git make
