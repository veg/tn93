
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



    min_overlap = atoi (argv[4]);
    if (min_overlap < 1) {
        cerr << "Minimum overlap must be a positive integer" << endl;
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
    
    while (fasta_result == 2) {
        fasta_result = readFASTA (F, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength, true);
        if (fasta_result == 1) {
            return 1;
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

