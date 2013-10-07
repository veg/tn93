
#include "tn93_shared.h"
#include <float.h>

#include "argparse_merge.hpp"

using namespace std;
using namespace argparse;


//---------------------------------------------------------------

void report_clusters (StringBuffer& defined_clusters, Vector& cluster_lengths, Vector& cluster_sizes, const long firstSequenceLength, const long minimum_cluster_size, bool json, char sep, FILE * out) {
    if (json) {
        bool do_comma = false;
        fprintf (out, "\n{");
        for (unsigned long k = 0; k < cluster_sizes.length(); k++) {
            if (cluster_sizes.value(k) >= minimum_cluster_size) {
                if (do_comma) {
                    fprintf (out, ",\n");
                } else {
                    fprintf (out, "\n");                
                }
                fprintf (out,"{\"read_count\":%ld, \"sequence\":\"", cluster_sizes.value(k));
                dump_fasta(stringText(defined_clusters, cluster_lengths, k), firstSequenceLength, out, false);
                fprintf (out,"\"}");
                do_comma = true;
            }
        } 
        fprintf (out, "\n}\n");    
    } else {
        for (unsigned long k = 0; k < cluster_sizes.length(); k++) {
            if (cluster_sizes.value(k) >= minimum_cluster_size) {
                fprintf (out, ">cluster_%ld%c%ld\n", k, sep, cluster_sizes.value(k));
                dump_fasta(stringText(defined_clusters, cluster_lengths, k), firstSequenceLength, out);
            }
        }
    }
} 

//---------------------------------------------------------------

void merge_current_clusters (StringBuffer& current_clusters, Vector& sequence_lengths, Vector& cluster_members, const long firstSequenceLength, const long min_overlap) {
    Vector cluster_merge;
    long   cluster_count = cluster_members.length();
    
    for (long k = 0; k < cluster_count; k++) {
        cluster_merge.appendValue(-1);
    }
    
    bool did_some_merges      = false;
    long merged_cluster_count = 0; 
    do {
       did_some_merges = false;
       for (long cluster_id = 1; cluster_id < cluster_count; cluster_id ++) {
            if (cluster_merge.value(cluster_id) >= 0) continue;
            char * cluster1 = stringText(current_clusters, sequence_lengths, cluster_id);
            long try_cluster = -1;
            
            #ifdef _OPENMP
              
            #if _OPENMP >= 200805 
              #pragma omp parallel shared(current_clusters, sequence_lengths, cluster_id, try_cluster, cluster_merge, cluster1, did_some_merges, merged_cluster_count)            
            #else 
              #pragma omp parallel shared(cluster_id, try_cluster, cluster_merge, cluster1, did_some_merges, merged_cluster_count)
            #endif 
            #endif
            {
                #pragma omp for schedule (dynamic)
                for (long cluster_id2 = 0; cluster_id2 < cluster_id; cluster_id2 ++) {
                    #pragma omp flush (try_cluster)
                    if (try_cluster < 0 && cluster_merge.value(cluster_id2) >= 0) {
                        char * cluster2 = stringText(current_clusters, sequence_lengths, cluster_id2);
                        if (perfect_match (cluster1, cluster2, firstSequenceLength) >= min_overlap) {
                          #pragma omp critical
                          {   
                                merge_two_sequences(cluster1, cluster2, firstSequenceLength);
                                did_some_merges = true;
                                cluster_merge.storeValue(cluster_id2, cluster_id);
                                merged_cluster_count ++;
                                try_cluster = cluster_id;
                           }
                           #pragma omp flush (try_cluster)
                        }   
                    }
                }
              }
       }
    
    } while (did_some_merges);
    
    if (merged_cluster_count) {
        StringBuffer compressed_clusters;
        Vector       compressed_lengths;
        Vector       compressed_members;
        compressed_lengths.appendValue (0);
        for (long k = 0; k < cluster_count; k++) {
            long merge_with = cluster_merge.value(k);
            if (merge_with < 0) {
                compressed_clusters.appendBuffer (stringText(current_clusters, sequence_lengths, k), firstSequenceLength);            
                compressed_clusters.appendChar (0);
                compressed_lengths.appendValue (compressed_clusters.length());
                compressed_members.appendValue(cluster_members.value (k));
            } else {
                compressed_members.storeValue (compressed_members.value(k) + compressed_members.value (merge_with), merge_with);
            }
        }
        compressed_clusters.swap (current_clusters);
        compressed_lengths.swap(sequence_lengths);
        compressed_members.swap (cluster_members);
        
        
    }
    
}

//---------------------------------------------------------------

void handle_a_sequence (StringBuffer& current_sequence, StringBuffer& current_clusters, Vector& sequence_lengths, Vector& cluster_members, const long firstSequenceLength, const long min_overlap, Vector* counts) {
    if (sequence_lengths.length() == 0) {
        sequence_lengths.appendValue(0);
    }

    unsigned long currently_defined_clusters = sequence_lengths.length()-1;
    long try_cluster = -1;
    
        #ifdef _OPENMP
          #if _OPENMP >= 200805 
            #pragma omp parallel for default(none) shared(currently_defined_clusters, try_cluster, sequence_lengths, current_sequence, current_clusters)
          #else 
            #pragma omp parallel for default(none) shared(currently_defined_clusters, try_cluster)
          #endif 
        #endif
        for (long cluster_index = 0; cluster_index < currently_defined_clusters; cluster_index ++) {
            #pragma omp flush (try_cluster)
            if (try_cluster < 0) {
                if (perfect_match (current_sequence.getString(), stringText(current_clusters, sequence_lengths, cluster_index), firstSequenceLength) >= min_overlap) {
                    #pragma omp critical
                    try_cluster = cluster_index;
                    #pragma omp flush (try_cluster)
                }
            }
        }
    
    if (try_cluster < 0) {
        current_clusters.appendBuffer(current_sequence.getString(), firstSequenceLength);
        current_clusters.appendChar (0);
        sequence_lengths.appendValue (current_clusters.length());
        cluster_members.appendValue(counts->value(0));
        
        if (cluster_members.length() % 1000 == 0) {
            merge_current_clusters (current_clusters, sequence_lengths, cluster_members, firstSequenceLength, min_overlap);
        }
        
    } else {
        cluster_members.storeValue (cluster_members.value(try_cluster) + counts->value(0),try_cluster);
        merge_two_sequences(current_sequence.getString(), stringText(current_clusters, sequence_lengths, try_cluster), firstSequenceLength);
    }
}   



//---------------------------------------------------------------

int main (int argc, const char * argv[])
{

    args_t args = args_t (argc, argv);
    initAlphabets();
    
    
    StringBuffer defined_clusters;
    
    Vector       cluster_lengths,
                 cluster_sizes,
                 copy_count;
    
    copy_count.appendValue(0);


    int resolutionOption;
    switch (args.ambig) {
      case resolve:
        resolutionOption = RESOLVE;
        break;
      case average:
        resolutionOption = AVERAGE;
        break;
      case skip:
        resolutionOption = SKIP;
        break;
      case gapmm:
        resolutionOption = GAPMM;
        break;
        
    }



    char automatonState = 0,
         fasta_result = 2;
    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence
 
    long sequences_read = 0,
         firstSequenceLength = 0;
          
    while (fasta_result == 2) {
        fasta_result = readFASTA (args.input, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength, true, &copy_count, args.counts_in_name);
        if (fasta_result == 1) {
            return 1;
        }
        handle_a_sequence (sequences, defined_clusters, cluster_lengths, cluster_sizes, firstSequenceLength, args.overlap, &copy_count);
        sequences_read ++;
        if (args.quiet == false && sequences_read % 1000 == 0) {
            cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << sequences_read << " reads (" << setw(8) << cluster_sizes.length() << " read-groups found)";
        }

    }
    merge_current_clusters (defined_clusters, cluster_lengths, cluster_sizes, firstSequenceLength,  args.overlap);
    
    if (args.quiet == false) {
      cerr << endl << "Found " << cluster_sizes.length() << " read-groups on " << sequences_read << " sequences." << endl;
    }
    
    
    /*
    long checksum = 0;
    for (long k = 0; k < cluster_sizes.length(); k++) {
        checksum += cluster_sizes.value(k);
    }
    cerr << sequences_read << ":" << checksum << endl;
    */

    report_clusters (defined_clusters, cluster_lengths, cluster_sizes, firstSequenceLength, args.cluster_size, args.json, args.counts_in_name, args.output);
    

    return 0;

}

