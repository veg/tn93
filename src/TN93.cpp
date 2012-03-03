
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string.h>

#include "stringBuffer.h"

#define  RESOLVE_A 0x01
#define  RESOLVE_C 0x02
#define  RESOLVE_G 0x04
#define  RESOLVE_T 0x08

using namespace std;

static char Usage[] = "TN93dist <FASTA file> <output file OR COUNT> <distance thershold> < how to handle ambiguities; one of RESOLVE, AVERAGE, SKIP> <output format; one of CSV, HYPHY> <minimum overlap between sequences: integer >= 1> [BOOTSTRAP 0 or 1] [SECOND FILE]",
ValidChars [] = "ACGTURYSWKMBDHVN?-",
empty      [] = "";

#define RESOLVE 0
#define AVERAGE 1
#define SKIP    2

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

#define GAP 17

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

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* prototypes */


/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}


/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}


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
        
        if (c1 == GAP || c2 == GAP)
            continue;
                
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
    
 	/*char	 nucs[]   = "ACGT",
    spacer[]  = "               ",
    spacer2[] = "----------------";
    
    cout << endl << totalNonGap << endl;
    
    fprintf (stderr,"\n\nPairiwse character counts\n*-%s|%s|%s|%s*\n", spacer2, spacer2, spacer2, spacer2);
    fprintf (stderr,"  %sA|%sC|%sG|%sT|\n",spacer,spacer,spacer,spacer); 
    fprintf (stderr,"*-%s|%s|%s|%s*\n", spacer2, spacer2, spacer2, spacer2);
    for (aux1 = 0; aux1 < 4; aux1++)
    {
        fprintf (stderr,"%c|",nucs[aux1]);
        for (aux2 = 0; aux2 < 4; aux2++)
        {
            fprintf (stderr,"%8g %6.2f%%|", pairwiseCounts[aux1][aux2], 100.*pairwiseCounts[aux1][aux2]/totalNonGap);
        }
        fprintf (stderr,"\n");
    }
    fprintf (stderr,"*-%s|%s|%s|%s*\n* ", spacer2, spacer2, spacer2, spacer2);
    for (aux1 = 0; aux1 < 4; aux1++)
    {
        fprintf (stderr,"%8g %6.2f%%|", nucFreq[aux1], 50.*nucFreq[aux1]/totalNonGap);
    }
    fprintf (stderr,"\n");
    fprintf (stderr,"*-%s|%s|%s|%s*\n* ", spacer2, spacer2, spacer2, spacer2);
    
    aux1 = totalNonGap-pairwiseCounts[0][0]-pairwiseCounts[1][1]-pairwiseCounts[2][2]-pairwiseCounts[3][3];
    fprintf (stderr, "\n\nHamming distance: %ld\np-distance: %g\n\n", aux1, aux1/totalNonGap);
	*/
    
	/*for (aux1 = 0; aux1 < 4; aux1++)
		for (aux2 = 0; aux2 < 4; aux2++)
		{
			nucFreq[aux1] += pairwiseCounts[aux1][aux2];
			nucFreq[aux2] += pairwiseCounts[aux1][aux2];
		}*/
	
	
    
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

int main (int argc, const char * argv[])
{
    if (argc != 7 && argc != 8 && argc != 9)
    {
        cerr << "Usage is `" << Usage << "'." << endl;
        exit(1);
    }
    
    const    char *S = argv[1];
    long     firstSequenceLength = 0,
    		 min_overlap = 1;
    		 
    double   distanceThreshold = atof (argv[3]);
    bool	 count_only 	   = strcmp ("COUNT", argv[2]) == 0;
    
    
    
    if (distanceThreshold < 0.0 || distanceThreshold > 1.0)
    {
        cerr << "Genetic distance threshold must be in [0,1]" << endl;
        return 1;
    }
    
    /*if (distanceThreshold > 0.)
    	min_overlap = 1./distanceThreshold;*/
    
    min_overlap = atoi (argv[6]);
    if (min_overlap < 1) {
        cerr << "Minimum overlap must be a positive integer" << endl;
        return 1;
    }
    
    FILE *F  = fopen(S, "r"),
         *FO = count_only ? NULL : fopen (argv[2], "w"),
         *F2 = argc == 9 ? fopen (argv[8], "r"): NULL;
    
    if (F == NULL)
    {
        cerr << "Cannot open file `" << S << "'." << endl;
        return 1;
    }
 
    if (FO == NULL && !count_only)
    {
        cerr << "Cannot open file `" << FO << "'." << endl;
        return 1;
    }
    
    if (F2 == NULL && argc == 9) {
        cerr << "Cannot open file `" << argv[8] << "'." << endl;
        return 1;
    
    }
    
    char resolutionOption = RESOLVE;
    
    if (strcmp (argv[4], "RESOLVE") == 0){
        resolutionOption = RESOLVE;
    } else if (strcmp (argv[4], "AVERAGE") == 0) {
        resolutionOption = AVERAGE;
    } else if (strcmp (argv[4], "SKIP") == 0) {
         resolutionOption = SKIP;
    }
 

    StringBuffer names,
    sequences;
    
    Vector       nameLengths,
    seqLengths;
    
    
    
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
    
    unsigned long seqLengthInFile1 = seqLengths.length()-1,
                  seqLengthInFile2 = 0;
    
    if (F2) {
        if (readFASTA (F2, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength))
            return 1;
        seqLengthInFile2 = seqLengths.length()-1-seqLengthInFile1;
        fclose (F2);
    }
    
    unsigned long sequenceCount = seqLengths.length()-1,
    pairwise      = argc == 9? seqLengthInFile2*seqLengthInFile1 : (sequenceCount-1) * (sequenceCount) / 2;
    
    if (argc < 9)
       cerr << "Read " << sequenceCount << " sequences of length " << firstSequenceLength << endl << "Will perform " << pairwise << " pairwise distance calculations";    
    else
        cerr << "Read " << seqLengthInFile1 << " sequences from file 1 and " << seqLengthInFile2 << " sequences from file 2 of length " << firstSequenceLength << endl << "Will perform " << pairwise << " pairwise distance calculations";
    
    double percentDone = 0.,
           max = 0.0;
    
    
    bool    doCSV = true;
    
    if (strcmp (argv[5], "HYPHY") == 0)
        doCSV = false;
    
    long * randFlag = NULL;
    if (argc == 8)
        if (atoi (argv[7]) > 0) {
            randFlag = new long [firstSequenceLength];
            init_genrand (time(NULL));
            for (long i = 0; i < firstSequenceLength; i++)
                randFlag[i] = genrand_int32 () % firstSequenceLength;
        }

       
    
    cerr << endl << "Progress: " << setw(8) << 0.0 << "% (" << setw(8) << 0 << " links found)";
    
    long pairIndex  = 0,
         foundLinks = 0;
         
    double *distanceMatrix = NULL;
    
    if (doCSV) {
        if (!count_only) {
            fprintf (FO, "Seq1,Seq2,Distance\n");
        }
    }
    else {
        distanceMatrix = new double [sequenceCount*sequenceCount];
        for (unsigned long i = 0; i < sequenceCount*sequenceCount; i++)
            distanceMatrix[i] = 100.;
        for (unsigned long i = 0; i < sequenceCount; i++)
            distanceMatrix[i*sequenceCount+i] = 0.;
    }
    
    long upperBound = argc == 9 ? seqLengthInFile1 : sequenceCount,
    	 skipped_comparisons = 0;
        
 #pragma omp parallel for default(none) shared(count_only, skipped_comparisons, resolutionOption, foundLinks,pairIndex,sequences,seqLengths,sequenceCount,firstSequenceLength,distanceThreshold, nameLengths, names, pairwise, percentDone,FO,cerr,max,randFlag,doCSV,distanceMatrix, upperBound, argc, seqLengthInFile1, seqLengthInFile2, min_overlap) 
    for (long seq1 = 0; seq1 < upperBound; seq1 ++)
    {
        char *n1 = stringText (names, nameLengths, seq1),
             *s1 = stringText(sequences, seqLengths, seq1);
        
        long lowerBound = argc == 9 ? seqLengthInFile1 : seq1 +1,  
             compsSkipped = 0;
        
        for (unsigned long seq2 = lowerBound; seq2 < sequenceCount; seq2 ++)
        {
            double thisD = computeTN93(s1, stringText(sequences, seqLengths, seq2), firstSequenceLength, resolutionOption, randFlag, min_overlap);
            
            if (thisD >= -1e-10 && thisD <= distanceThreshold)
            {
                #pragma omp critical
                foundLinks ++;
                //char *s2 = stringText(sequences, seqLengths, seq1);
                if (!count_only){
					if (doCSV){
						#pragma omp critical
						fprintf (FO,"%s,%s,%g\n", n1, stringText (names, nameLengths, seq2), thisD);
					} else {
						distanceMatrix[seq1*sequenceCount+seq2] = thisD;
						distanceMatrix[seq2*sequenceCount+seq1] = thisD;
					}
				}
            }
            if (thisD > max)
            {
                #pragma omp critical
                max = thisD;
            }
            else
            {
            	if (thisD <= -0.5) {
                	
            		compsSkipped ++;
            	}
            }
            
            
        }
        #pragma omp critical
        {
        skipped_comparisons += compsSkipped;
        pairIndex += (argc < 9) ? (sequenceCount - seq1 - 1) : seqLengthInFile2;
        }
        
        if (pairIndex * 100. / pairwise - percentDone > 0.1 || seq1 == (long)sequenceCount - 1)
        {
            #pragma omp critical
            percentDone = pairIndex * 100. / pairwise;
            #pragma omp critical
            cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << percentDone << "% (" << setw(8) << foundLinks << " links found)";
        }
   }

    if (!count_only && !doCSV)
    {
        fprintf (FO, "{");
        for (unsigned long seq1 = 0; seq1 < sequenceCount; seq1 ++)
        {
            fprintf (FO, "\n{%g", distanceMatrix[seq1*sequenceCount]);
            for (unsigned long seq2 = 1; seq2 < sequenceCount; seq2++)
            {
                fprintf (FO, ",%g", distanceMatrix[seq1*sequenceCount+seq2]);
            }
            fprintf (FO, "}");
        }
        fprintf (FO, "\n}");

        delete [] distanceMatrix;
    }
    
    cerr << endl;        
    cerr << "Actual comparisons performed " << pairwise-skipped_comparisons << endl;
    cerr << "Maximum distance = " << max << endl;
    
    if (count_only) {
    	cout << "Found " << foundLinks << " links among " << pairwise-skipped_comparisons << " pairwise comparisons" << endl;
    }
    
    if (randFlag)
        delete [] randFlag;
    
    if (FO)
    	fclose (FO);
    return 0;
    
}

