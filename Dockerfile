# Docker image for a tn93 development environment
FROM oraclelinux:8

# Set up environment and install dependencies
RUN yum -y update && \
    yum install -y cmake gcc-c++ git make

    # The current tn93 OpenMP syntax is incompatible with GCC 8.5.0, which is what gets installed by the above commands.
    # To compile tn93, make the following two edits to the tn93 source code:
    #sed -i 's/shared(currently_defined_clusters, try_cluster, sequence_lengths, current_sequence, current_clusters, firstSequenceLength, min_overlap)/shared(currently_defined_clusters, try_cluster, sequence_lengths, current_sequence, current_clusters)/g' src/read_reducer.cpp && \
    #sed -i 's/shared(my_distance_estimate,nodeParents,workingNodes,distanceEstimates, step_penalty, min_overlap, resolutionOption, firstSequenceLength, theSequence, left_to_do)/shared(my_distance_estimate,nodeParents,workingNodes,distanceEstimates)/g' src/ShortestPathTN93.cpp
