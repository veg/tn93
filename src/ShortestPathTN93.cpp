
#include "tn93_shared.h"
#include <float.h>

using namespace std;

static char Usage[] = "ShortestPathTN93"
                      "\n\t<FASTA file OR - for stdin >"
                      "\n\t<output file OR - for stdout>"
                      "\n\t< how to handle ambiguities; one of RESOLVE, AVERAGE, SKIP, GAPMM>"
                      "\n\t<minimum overlap between sequences: integer >= 1>"
                      "\n\t<0-based index of the source sequence"
                      "\n\t<output format: FASTA or JSON>"
                      "\n\t<large step penalty: positive real number>"
                      "\n\t<report estimate paths for these sequences -- 0-based indices>";



//---------------------------------------------------------------

double    computeTransformedTN93 (const unsigned long seq1, const unsigned long seq2, const long firstSequenceLength, const char resolutionOption, const long min_overlap, const double step_penalty, double *ds = NULL) {

    char *s1 = stringText (sequences, seqLengths, seq1);

    double thisD = computeTN93(s1, stringText(sequences, seqLengths, seq2), firstSequenceLength, resolutionOption, NULL, min_overlap),
           d = thisD < 0.?1.e100:exp(thisD*step_penalty)-1.;


    if (ds)
        *ds = thisD;

    return d;
}



//---------------------------------------------------------------

void initializeSingleSource (const unsigned long seq_count, const unsigned long source) {
    for (unsigned long idx = 0; idx < seq_count; idx++) {
        distanceEstimates.appendValue (DBL_MAX);
        workingNodes.appendValue (idx);
        nodeParents.appendValue  (-1);
    }
    distanceEstimates.storeValue (0., source);
}


//---------------------------------------------------------------

void reportPathToSource (const unsigned long which_index, FILE* output, const bool is_json, const long firstSequenceLength, const char resolutionOption, const long min_overlap, const double step_penalty) {
    if (which_index < nodeParents.length()) {
        if (is_json) {
            char * sname = stringText (names, nameLengths, which_index);
            fprintf (output, "\"%s\": [\"%s\"",  sname, sname);
        } else {
            dump_sequence_fasta (which_index, output, firstSequenceLength);
        }
        long    last_index    = which_index,
                current_index = nodeParents.value(which_index);

        while (current_index >= 0) {
            if (is_json) {
                fprintf (output, ",\"%s\"",  stringText (names, nameLengths, current_index));
            } else {
                double d [2];
                d[1] = computeTransformedTN93 (current_index, last_index, firstSequenceLength, resolutionOption, min_overlap, step_penalty, d);
                dump_sequence_fasta (current_index, output,firstSequenceLength, d);
            }
            last_index = current_index;
            current_index = nodeParents.value(current_index);
        }
        if (is_json) {
            fprintf (output, "]");
        } else {
            fprintf (output, "\n");
        }

    } else {
        cerr << "Sequence index " << which_index << " is out of range." << endl;
    }
}

//---------------------------------------------------------------

void relaxDistanceEstimates (const unsigned long theSequence, const long firstSequenceLength, const char resolutionOption, const long min_overlap, const double step_penalty) {
    const unsigned long left_to_do           = workingNodes.length();
    double               my_distance_estimate = distanceEstimates.value (theSequence);

    #pragma omp parallel for default(none) shared(my_distance_estimate,nodeParents,workingNodes,distanceEstimates)

    for (long remaining = 0; remaining < left_to_do; remaining ++) {
        const unsigned long working_index = workingNodes.value(remaining);
        double new_estimate = computeTransformedTN93 (theSequence, working_index, firstSequenceLength, resolutionOption, min_overlap, step_penalty);
        if (new_estimate < DBL_MAX) {
            if (new_estimate + my_distance_estimate < distanceEstimates.value (working_index)) {
                distanceEstimates.storeValue (new_estimate + my_distance_estimate, working_index);
                nodeParents.storeValue (theSequence, working_index);
            }
        }
    }
}


//---------------------------------------------------------------

int main (int argc, const char * argv[])
{
    if (argc < 9)
    {
        cerr << "Usage is `" << Usage << "'." << endl;
        exit(1);
    }

    const    char *S = argv[1];

    long     firstSequenceLength = 0,
             min_overlap = 1;


    /*if (distanceThreshold > 0.)
    	min_overlap = 1./distanceThreshold;*/

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




    char automatonState = 0;
    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence

    nameLengths.appendValue (0);
    seqLengths.appendValue (0);
    initAlphabets ();

    if (readFASTA (F, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength) == 1)
        return 1;

    fclose(F);

    unsigned long sequenceCount = seqLengths.length()-1;

    cerr << "Read " << sequenceCount << " sequences of length " << firstSequenceLength << endl ;

    unsigned long source = atoi (argv[5]);
    if (source < 0 || source >= sequenceCount) {
        cerr << "Invalid source sequence index";
        return 1;
    } else {
        cerr << "Using " << stringText (names, nameLengths, source) << " as the source" << endl;
    }

    initializeSingleSource (sequenceCount, source);

    double percentDone = 0.0,
           normalizer  = 100./sequenceCount,
           step_penalty = 10.;

    step_penalty = atof (argv[7]);
    if (step_penalty <= 0.0) {
        cerr << "Invalid step penalty";
        return 1;
    } else {
        cerr << "Using step penalty of " << step_penalty << endl;
    }

    while (workingNodes.length()) {
        unsigned long add_this_node = workingNodes.extractMin (distanceEstimates);
        relaxDistanceEstimates (add_this_node,firstSequenceLength,resolutionOption,min_overlap,step_penalty);
        if ((sequenceCount-workingNodes.length()) * normalizer - percentDone > 0.1 || workingNodes.length () == 0) {
            cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << percentDone << "% done ";
            percentDone = (sequenceCount-workingNodes.length()) * normalizer;
        }
    }
    cerr << endl;

    bool is_json = strcmp (argv[6],"JSON") == 0;

    if (is_json) {
        fprintf (FO, "\n{\n");
    }

    for (long which_arg = 8; which_arg < argc; which_arg ++) {
        reportPathToSource (atoi (argv[which_arg]),FO, is_json, firstSequenceLength,resolutionOption,min_overlap,step_penalty);
        if (is_json && which_arg < argc - 1) {
            fprintf (FO, ",\n");
        }
    }

    if (is_json) {
        fprintf (FO, "\n}\n");
    }

    if (FO != stdout)
        fclose (FO);
    return 0;

}

