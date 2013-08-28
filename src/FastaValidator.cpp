#include "tn93_shared.h"

using namespace std;

static char Usage[] = "Validate FASTA"
                      "\n\t<FASTA file OR - for stdin >";

//---------------------------------------------------------------

int main (int argc, const char * argv[])
{
    if (argc != 2)
    {
        cerr << "Usage is `" << Usage << "'." << endl;
        exit(1);
    }

    const    char *S = argv[1];
    long     firstSequenceLength = 0;

    FILE *F  = strcmp (S,"-") == 0 ? stdin: fopen(S, "r");

    if (F == NULL)
    {
        cerr << "Cannot open file `" << S << "'." << endl;
        return 1;
    }

    initAlphabets();
    if (validateFASTA (F) == 1)
        return 1;

    fclose(F);
    return 0;
}

