#include "tn93_shared.h"

using namespace std;

static char Usage[] = "Validate FASTA"
                      "\n\t<FASTA file OR - for stdin >";

//---------------------------------------------------------------

int main (int argc, const char * argv[])
{
    if (argc <= 1)
    {
        cerr << "Usage is `" << Usage << "'." << endl;
        exit(1);
    }

    const char *S = argv[1];
    long firstSequenceLength = 0;

    FILE *F = strcmp (S,"-") == 0 ? stdin: fopen(S, "r");

    if (F == NULL)
    {
        cerr << "Cannot open file `" << S << "'." << endl;
        return 1;
    }

    StringBuffer names, sequences;
    Vector nameLengths, seqLengths;
    char automatonState = 0;

    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence
    nameLengths.appendValue (0);
    seqLengths.appendValue (0);
    initAlphabets();
    if (readFASTA (F, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength) == 1)
        return 1;

    fclose(F);
    cout << "Valid FASTA file" << endl;
    return 0;
}
