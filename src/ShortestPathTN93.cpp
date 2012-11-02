
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <cfloat>

#include "stringBuffer.h"

#define  RESOLVE_A 0x01
#define  RESOLVE_C 0x02
#define  RESOLVE_G 0x04
#define  RESOLVE_T 0x08

using namespace std;

static char Usage[] = "ShortestPathTN93"
                      "\n\t<FASTA file OR - for stdin >"
                      "\n\t<output file OR - for stdout>"
                      "\n\t< how to handle ambiguities; one of RESOLVE, AVERAGE, SKIP, GAPMM>"
                      "\n\t<minimum overlap between sequences: integer >= 1>"
                      "\n\t<0-based index of the source sequence"
                      "\n\t<output format: FASTA or JSON>"
                      "\n\t<report estimate paths for these sequences -- 0-based indices>",


ValidChars [] = "ACGTURYSWKMBDHVN?-",
empty      [] = "";

#define RESOLVE 0
#define AVERAGE 1
#define SKIP    2
#define GAPMM   3

static unsigned char   resolutions [] = { RESOLVE_A, 
    RESOLVE_C,
    RESOLVE_G,
    RESOLVE_T,
    RESOLVE_T, // U - 4
    RESOLVE_A | RESOLVE_G, // R - 5
    RESOLVE_C | RESOLVE_T, // Y - 6
    RESOLVE_C | RESOLVE_G, // S - 7
    RESOLVE_A | RESOLVE_T, // W - 8
    RESOLVE_G | RESOLVE_T, // K - 9
    RESOLVE_A | RESOLVE_C, // M - 10
    RESOLVE_C | RESOLVE_G | RESOLVE_T, // B - 11
    RESOLVE_A | RESOLVE_G | RESOLVE_T, // D - 12
    RESOLVE_A | RESOLVE_C | RESOLVE_T, // H - 13
    RESOLVE_A | RESOLVE_C | RESOLVE_G, // V - 14
    RESOLVE_A | RESOLVE_C | RESOLVE_G | RESOLVE_T , // N - 15
    RESOLVE_A | RESOLVE_C | RESOLVE_G | RESOLVE_T , // ? - 16
    0. };

#define N_CHAR 15
#define GAP    17

static  double   resolutionsCount [] = { 1.f,
    1.f,
    1.f,
    1.f,
    1.f,
    1./2.f, // R
    1./2.f, // Y
    1./2.f, // S
    1./2.f, // S
    1./2.f, // W
    1./2.f, // K
    1./2.f, // M
    1./3.f, // B
    1./3.f, // D
    1./3.f, // H
    1./3.f, // V
    1./4.f, // N
    1./4.f, // ?
    0.f };

char validFlags [256];

//---------------------------------------------------------------

StringBuffer names,
             sequences;

Vector       nameLengths,
             seqLengths,
             workingNodes,
             nodeParents;

VectorDouble distanceEstimates;

long         firstSequenceLength = 0,
             min_overlap = 1;

char         resolutionOption = RESOLVE;


//---------------------------------------------------------------

long stringLength (Vector& lengths, unsigned long index)
{
    if (index < lengths.length() - 1)
        return lengths.value(index+1) - lengths.value(index) - 1;
    
    return -1;
}

//---------------------------------------------------------------

char* stringText (StringBuffer& strings, Vector& lengths, unsigned long index)
{
    if (index < lengths.length() - 1)
        return strings.getString() + lengths.value(index);
    
    return empty;
}

/*---------------------------------------------------------------------------------------------------- */

double		computeTN93 (char * s1, char *s2,  unsigned long L, char matchMode, long* randomize, long min_overlap)
{
	char useK2P   = 0;
    
	long aux1;
    
	double auxd,
    nucFreq[4]		= {0,0,0,0},
    fY,
    fR,
    K1,
    K2,
    K3,
    AG,
    CT,
    ti,
    tv,
    totalNonGap	= 0.,
    nucF[4],
    pairwiseCounts [4][4];
    
    for (long i = 0; i < 4; i++)
        for (long j = 0; j < 4; j++)
            pairwiseCounts[i][j] = 0.;
            
    
    for (unsigned long p = 0; p < L; p++)
    {
        unsigned char c1, c2;
        
        if (randomize) {
            long pi = randomize[p];
            c1 = s1[pi];
            c2 = s2[pi];
        } else {
            c1 = s1[p];
            c2 = s2[p];
        }
        
        if (c1 == GAP || c2 == GAP) {  
            if (matchMode == GAPMM) {
                if (c1 == GAP && c2 == GAP)
                   continue;
                else {
                    if (c1 == GAP) {
                        c1 = N_CHAR;
                    } else {
                        c2 = N_CHAR;
                    }
                }
            } else
                continue;
        }
                
        if (c1 < 4)
        {
            if (c2 < 4)
            {
                pairwiseCounts [c1][c2] += 1.; 
            }
            else
            {
                if (matchMode != SKIP)
                {
                    if (resolutionsCount[c2] > 0.)
                    {
                        if (matchMode == RESOLVE)
                            if (resolutions[c2] & (1 << c1))
                            {
                                pairwiseCounts[c1][c1] += 1.;
                                continue;
                            }
                        if (resolutions[c2] & RESOLVE_A)
                            pairwiseCounts[c1][0] += resolutionsCount[c2]; 
                        if (resolutions[c2] & RESOLVE_C)
                            pairwiseCounts[c1][1] += resolutionsCount[c2]; 
                        if (resolutions[c2] & RESOLVE_G)
                            pairwiseCounts[c1][2] += resolutionsCount[c2]; 
                        if (resolutions[c2] & RESOLVE_T)
                            pairwiseCounts[c1][3] += resolutionsCount[c2]; 
                    }
                }
            }
        }
        else
        {
            if (matchMode != SKIP)
            {
                if (c2 < 4)
                {
                    if (resolutionsCount[c1] > 0.)
                    {
                        if (matchMode == RESOLVE)
                            if (resolutions[c1] & (1 << c2))
                            {
                                pairwiseCounts[c2][c2] += 1.;
                                continue;
                            }
                        
                        if (resolutions[c1] & RESOLVE_A)
                            pairwiseCounts[0][c2] += resolutionsCount[c1]; 
                        if (resolutions[c1] & RESOLVE_C)
                            pairwiseCounts[1][c2] += resolutionsCount[c1]; 
                        if (resolutions[c1] & RESOLVE_G)
                            pairwiseCounts[2][c2] += resolutionsCount[c1]; 
                        if (resolutions[c1] & RESOLVE_T)
                            pairwiseCounts[3][c2] += resolutionsCount[c1]; 
                    }            
                }
                else // ambig and ambig
                {
                    double norm = resolutionsCount[c1] * resolutionsCount[c2];
                    //cout << int(c1) << ":" << int(c2) << "/" << norm << endl;
                    if (norm > 0.0)
                    {
                        if (matchMode == RESOLVE)
                        {
                            
                        }
                        
                        char indexer = 1;
                        for (long i = 0; i < 4; i ++, indexer <<= 1)
                        {
                            if (resolutions[c1] & indexer)
                            {
                                char indexer2 = 1;
                                for (long j = 0; j < 4; j ++, indexer2 <<= 1)
                                {
                                    if (resolutions [c2] & indexer2)
                                    {
                                        pairwiseCounts[i][j] += norm;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
     
    for (long c1 = 0; c1 < 4; c1++)
        for (long c2 = 0; c2 < 4; c2++)
        {
            totalNonGap   += pairwiseCounts[c1][c2];
            nucFreq [c1]  += pairwiseCounts[c1][c2];
            nucFreq [c2]  += pairwiseCounts[c1][c2];
        }
    
    if (totalNonGap <= min_overlap) {
    	return -1.;
    }
       
	totalNonGap = 2./(nucFreq[0] + nucFreq[1] + nucFreq[2] + nucFreq[3]);
       
	auxd = 1./(nucFreq[0] + nucFreq[1] + nucFreq[2] + nucFreq[3]);
	for (aux1 = 0; aux1 < 4; aux1++)
		nucF[aux1] = nucFreq[aux1]*auxd;
	
	fR = nucF[0]+nucF[2];
	fY = nucF[1]+nucF[3];
	
	if (nucFreq[0] == 0 || nucFreq[1] == 0 || nucFreq[2] == 0 || nucFreq[3] == 0)
		useK2P = 1;
	else
	{
		K1 = 2.*nucF[0]*nucF[2]/fR;
		K2 = 2.*nucF[1]*nucF[3]/fY;
		K3 = 2.*(fR*fY - nucF[0]*nucF[2]*fY/fR - nucF[1]*nucF[3]*fR/fY);
	}
    
	AG = (pairwiseCounts[0][2] + pairwiseCounts[2][0])*totalNonGap;
	CT = (pairwiseCounts[1][3] + pairwiseCounts[3][1])*totalNonGap;
	tv = 1.-((pairwiseCounts[0][0] + pairwiseCounts[1][1] + pairwiseCounts[2][2] + pairwiseCounts[3][3])*totalNonGap + 
			 AG+CT);
    
    
	if (useK2P)
	{
	    ti = AG+CT;
		AG = 1.-2.*ti-tv;
		CT = 1.-2.*tv;
		if (AG > 0. && CT > 0.)
			return -0.5*log(AG)-0.25*log(CT);
	}
	else
	{
		AG	= 1.-AG/K1 - 0.5*tv/fR;
		CT	= 1.-CT/K2 - 0.5*tv/fY;
		tv  = 1.-0.5 * tv/fY/fR;
		if (AG > 0. && CT > 0. && tv > 0)
			return - K1*log(AG) - K2*log(CT) - K3 * log (tv);
	}
	
	return 1000.0;
}


//---------------------------------------------------------------

void addASequenceToList (StringBuffer& sequences, Vector& seqLengths, long &firstSequenceLength, StringBuffer& names, Vector& nameLengths)
{
    sequences.appendChar ('\0');
    seqLengths.appendValue (sequences.length());
    if (seqLengths.length() == 2)
    {
        firstSequenceLength = stringLength (seqLengths, 0);
        
        if (firstSequenceLength <= 0)
        {
            cerr << "First sequence length must be positive." << endl;
            exit (1);
        }
    }
    else
    {
        if (stringLength (seqLengths, seqLengths.length()-2) != firstSequenceLength)
        {
            cerr << "All sequences must have the same length (" << firstSequenceLength << "), but sequence '" << stringText (names, nameLengths, nameLengths.length()-2) << "' had length " << stringLength (seqLengths, seqLengths.length()-2);
            exit (1);                           
        }
    }
}

//---------------------------------------------------------------

int readFASTA (FILE* F, char& automatonState,  StringBuffer &names, StringBuffer& sequences, Vector &nameLengths, Vector &seqLengths, long& firstSequenceLength)
{
    while (1)
    {
   		int currentC = fgetc (F);	
        if (feof (F))
            break;
        switch (automatonState)
        {
            case 0:
            {
   				if (currentC == '>' || currentC == '#')
   					automatonState = 1;
   				break;
            }
            case 1:
            {
   				if (currentC == '\n' || currentC == '\r')
                {
                    names.appendChar   ('\0');
                    nameLengths.appendValue (names.length());
                    //cout << stringText (names, nameLengths, nameLengths.length () - 2) << endl;
                    if (stringLength (nameLengths, nameLengths.length()-2) <= 0)
                    {
                        cerr << "Sequence names must be non-empty." << endl;
                        return 1;
                    }
                    
                    automatonState = 2;
                }
                else
                {
                    names.appendChar (currentC);
                }
   				break;
            }
            case 2:
            {
                currentC = toupper (currentC);
   				if (validFlags [currentC] >= 0)
                {
                    //cout << "Append " << currentC << endl;
                    sequences.appendChar (validFlags [currentC]);
                }
   				else
   					if (currentC == '>' || currentC == '#')
                    {
                        addASequenceToList (sequences, seqLengths, firstSequenceLength, names, nameLengths);
                        automatonState = 1;
                    }
                break;
            }
        }
        
    }
    
    if (automatonState == 2){
        addASequenceToList (sequences, seqLengths, firstSequenceLength, names, nameLengths);
        automatonState = 1;
    } else {
        cerr << "Unexpected end of file" << endl;
        return 1;
    }
    return 0;
}



//---------------------------------------------------------------

double    computeTransformedTN93 (unsigned long seq1, unsigned long seq2) {

    char *n1 = stringText (names, nameLengths, seq1),
         *s1 = stringText (sequences, seqLengths, seq1);
    
    double thisD = computeTN93(s1, stringText(sequences, seqLengths, seq2), firstSequenceLength, resolutionOption, NULL, min_overlap),
    d = exp(thisD*20.)-1.;
    
    
    //#pragma omp critical
    //cout << thisD << " -> " << d << endl;
    return d;
    //return (long)(exp(10.*thisD)-1.0)*10000;
}



//---------------------------------------------------------------

void initializeSingleSource (unsigned long seq_count, unsigned long source) {
    for (unsigned long idx = 0; idx < seq_count; idx++) {
        distanceEstimates.appendValue (DBL_MAX);
        workingNodes.appendValue (idx);
        nodeParents.appendValue  (-1);
    }
    distanceEstimates.storeValue (0., source);
}

//---------------------------------------------------------------

void dump_sequence_fasta (unsigned long index, FILE* output) {
    fprintf (output, ">%s\n", stringText (names, nameLengths, index));
    char *s1 = stringText (sequences, seqLengths, index);
    for (long c = 0; c < firstSequenceLength; c++) {
        fputc( ValidChars[s1[c]], output);
    }
    fprintf (output, "\n");
}

//---------------------------------------------------------------

void reportPathToSource (const unsigned long which_index, FILE* output, bool is_json) {
    if (which_index < nodeParents.length()) {
        if (is_json) {
            char * sname = stringText (names, nameLengths, which_index);
            fprintf (output, "\"%s\": [\"%s\"",  sname, sname);
        } else {
            dump_sequence_fasta (which_index, output);
        }
        long current_index = nodeParents.value(which_index);
        
        while (current_index >= 0) {
            if (is_json) {
                fprintf (output, ",\"%s\"",  stringText (names, nameLengths, current_index));
            } else {
                dump_sequence_fasta (current_index, output);
            }
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

void relaxDistanceEstimates (unsigned long theSequence) {
    const unsigned long left_to_do           = workingNodes.length();
    double               my_distance_estimate = distanceEstimates.value (theSequence);
                        
 #pragma omp parallel for default(none) shared(my_distance_estimate,nodeParents,workingNodes,theSequence,distanceEstimates ) 
 
    for (long remaining = 0; remaining < left_to_do; remaining ++) {
        const unsigned long working_index = workingNodes.value(remaining);
        double new_estimate = computeTransformedTN93 (theSequence, working_index);
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
    if (argc < 8)
    {
        cerr << "Usage is `" << Usage << "'." << endl;
        exit(1);
    }
    
    const    char *S = argv[1];
    		 
    
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
    
    
    if (strcmp (argv[3], "GAPMM") == 0){
        resolutionOption = GAPMM;
    } else if (strcmp (argv[3], "AVERAGE") == 0) {
        resolutionOption = AVERAGE;
    } else if (strcmp (argv[3], "SKIP") == 0) {
         resolutionOption = SKIP;
    }
 
    
    
    
    for (int i = 0; i < 256; i++)
        validFlags [i] = -1;
    
    for (unsigned int i = 0; i < strlen (ValidChars); i++)
        validFlags [(unsigned char)ValidChars[i]] = i;
    
    char automatonState = 0;
    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence
    
    nameLengths.appendValue (0);
    seqLengths.appendValue (0);
    
    if (readFASTA (F, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength))
        return 1;
        
    fclose(F);
        
    unsigned long sequenceCount = seqLengths.length()-1;
    
    cerr << "Read " << sequenceCount << " sequences of length " << firstSequenceLength << endl ;
    
    unsigned long source = atoi (argv[5]);
    if (source < 0 || source >= sequenceCount) {
        cerr << "Invalid source sequence index";
        return 1;
    }
    
    initializeSingleSource (sequenceCount, source);
    
    double percentDone = 0.0,
           normalizer  = 100./sequenceCount;
    
    while (workingNodes.length()) {
        unsigned long add_this_node = workingNodes.extractMin (distanceEstimates);
        relaxDistanceEstimates (add_this_node);
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
    
    for (long which_arg = 7; which_arg < argc; which_arg ++) {
        reportPathToSource (atoi (argv[which_arg]),FO, is_json);
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

