# Docker image for a tn93 development environment (Red Hat Universal Base Image 8.10)
FROM redhat/ubi8:8.10

# Set up environment and install dependencies
RUN yum -y update && \
    yum install -y cmake gcc-c++ gcc-toolset-12 git && \
    echo 'source /opt/rh/gcc-toolset-12/enable' > ~/.bashrc && \
    source ~/.bashrc

# To compile tn93 within the development environment:
#   cmake .
#   make
