# Docker image for tn93
FROM ubuntu:20.04

# Set up environment and install dependencies
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y cmake g++ git make
