
#include "tn93_shared.h"
#include <float.h>

using namespace std;

static char Usage[] = "read_reducer"
                      "\n\t<FASTA file OR - for stdin >"
                      "\n\t<output file OR - for stdout>"
                      "\n\t<how to handle ambiguities; one of RESOLVE, AVERAGE, SKIP, GAPMM>"
                      "\n\t<minimum overlap between sequences: integer >= 1>"
                      "\n\t<output format: FASTA or JSON>"
                      "\n\t<report only clusters with minimum of this many reads: integer >= 1>";


//---------------------------------------------------------------

void report_clusters (StringBuffer& defined_clusters, Vector& cluster_lengths, Vector& cluster_sizes, const long firstSequenceLength, const long minimum_cluster_size, bool json, FILE * out) {
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
                fprintf (out,"{read_count:%ld, sequence:\"", cluster_sizes.value(k));
                dump_fasta(stringText(defined_clusters, cluster_lengths, k), firstSequenceLength, out, false);
                fprintf (out,"\"}");
                do_comma = true;
            }
        } 
        fprintf (out, "\n}\n");    
    } else {
        for (unsigned long k = 0; k < cluster_sizes.length(); k++) {
            if (cluster_sizes.value(k) >= minimum_cluster_size) {
                fprintf (out, ">cluster_%ld_%ld\n", k, cluster_sizes.value(k));
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
            
            #pragma omp parallel for default(none) shared(sequence_lengths, current_clusters, cluster_id, try_cluster, cluster_merge, cluster1, did_some_merges, merged_cluster_count)
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
    /*
    for (long k = 0; k < cluster_count; k++) {
        printf ("%ld => %ld\n", k, cluster_merge.value(k));
    }
    */
    
}

//---------------------------------------------------------------

void handle_a_sequence (StringBuffer& current_sequence, StringBuffer& current_clusters, Vector& sequence_lengths, Vector& cluster_members, const long firstSequenceLength, const long min_overlap) {
    if (sequence_lengths.length() == 0) {
        sequence_lengths.appendValue(0);
    }

    unsigned long currently_defined_clusters = sequence_lengths.length()-1;
    long try_cluster = -1;
    
    
    #pragma omp parallel for default(none) shared(current_clusters, currently_defined_clusters, try_cluster, current_sequence, sequence_lengths)
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
        cluster_members.appendValue(1);
        
        if (cluster_members.length() % 1000 == 0) {
            merge_current_clusters (current_clusters, sequence_lengths, cluster_members, firstSequenceLength, min_overlap);
        }
        
    } else {
        cluster_members.storeValue (cluster_members.value(try_cluster) + 1,try_cluster);
        merge_two_sequences(current_sequence.getString(), stringText(current_clusters, sequence_lengths, try_cluster), firstSequenceLength);
    }
}   



//---------------------------------------------------------------

int main (int argc, const char * argv[])
{
    if (argc != 7)
    {
        cerr << "Usage is `" << Usage << "'." << endl;
        exit(1);
    }

    const    char *S = argv[1];

    long     firstSequenceLength  = 0,
             min_overlap          = 50,
             minimum_cluster_size = 5;
             
    StringBuffer defined_clusters;
    Vector       cluster_lengths,
                 cluster_sizes;


    min_overlap = atoi (argv[4]);
    if (min_overlap < 1) {
        cerr << "Minimum overlap must be a positive integer" << endl;
        return 1;
    }

    minimum_cluster_size = atoi (argv[6]);
    if (minimum_cluster_size < 1) {
        cerr << "Minimum cluster size must be a positive integer" << endl;
        return 1;
    }
    
    FILE *F  = strcmp (S,"-") == 0 ? stdin: fopen(S, "r"),
          *FO = strcmp (argv[2],"-") == 0? stdout: fopen (argv[2], "w");

    if (F == NULL) {
        cerr << "Cannot open file `" << S << "'." << endl;
        return 1;
    }

    if (FO == NULL) {
        cerr << "Cannot open file `" << FO << "'." << endl;
        return 1;
    }

    char resolutionOption = RESOLVE;

    if (strcmp (argv[3], "GAPMM") == 0) {
        resolutionOption = GAPMM;
    } else if (strcmp (argv[3], "AVERAGE") == 0) {
        resolutionOption = AVERAGE;
    } else if (strcmp (argv[3], "SKIP") == 0) {
        resolutionOption = SKIP;
    }




    char automatonState = 0,
         fasta_result = 2;
    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence

    initAlphabets();
 
    long sequences_read = 0;
          
    while (fasta_result == 2) {
        fasta_result = readFASTA (F, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength, true);
        if (fasta_result == 1) {
            return 1;
        }
        handle_a_sequence (sequences, defined_clusters, cluster_lengths, cluster_sizes, firstSequenceLength, min_overlap);
        sequences_read ++;
        if (sequences_read % 1000 == 0) {
            cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << sequences_read << " reads (" << setw(8) << cluster_sizes.length() << " read-groups found)";
        }

    }
    merge_current_clusters (defined_clusters, cluster_lengths, cluster_sizes, firstSequenceLength, min_overlap);
    cerr << endl << "Found " << cluster_sizes.length() << " read-groups on " << sequences_read << " sequences." << endl;
    
    fclose(F);
    
    long checksum = 0;
    for (long k = 0; k < cluster_sizes.length(); k++) {
        checksum += cluster_sizes.value(k);
    }
    cerr << sequences_read << ":" << checksum << endl;

    report_clusters (defined_clusters, cluster_lengths, cluster_sizes, firstSequenceLength, minimum_cluster_size, strcmp (argv[5],"JSON") == 0, FO);
    

    if (FO != stdout)
        fclose (FO);
    return 0;

}

