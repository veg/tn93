
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

void handle_a_sequence (StringBuffer& current_sequence, StringBuffer& current_clusters, Vector& sequence_lengths, Vector& cluster_members, const long firstSequenceLength, const long min_overlap) {
    if (sequence_lengths.length() == 0) {
        sequence_lengths.appendValue(0);
    }

    bool make_a_new_cluster = true;
    unsigned long currently_defined_clusters = sequence_lengths.length()-1,
                  try_cluster;
    
    
    for (try_cluster = 0; try_cluster < currently_defined_clusters; try_cluster ++) {
        if (perfect_match (current_sequence.getString(), stringText(current_clusters, sequence_lengths, try_cluster), firstSequenceLength) >= min_overlap) {
            make_a_new_cluster = false;
            break;
        }
    }
    
    if (make_a_new_cluster) {
        current_clusters.appendBuffer(current_sequence.getString(), firstSequenceLength);
        current_clusters.appendChar (0);
        sequence_lengths.appendValue (current_clusters.length());
        cluster_members.appendValue(1);
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
    cerr << endl;
    
    for (long k = 0; k < cluster_sizes.length(); k++) {
        if (cluster_sizes.value(k) >= minimum_cluster_size) {
            printf (">cluster%ld_%ld\n", k, cluster_sizes.value(k));
            dump_fasta(stringText(defined_clusters, cluster_lengths, k), firstSequenceLength, stdout);
        }
    }

    fclose(F);

    
    bool is_json = strcmp (argv[6],"JSON") == 0;

    if (is_json) {
        fprintf (FO, "\n{\n");
    }

    
    if (is_json) {
        fprintf (FO, "\n}\n");
    }

    if (FO != stdout)
        fclose (FO);
    return 0;

}

