
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "tn93_shared.h"

#define  RESOLVE_A 0x01
#define  RESOLVE_C 0x02
#define  RESOLVE_G 0x04
#define  RESOLVE_T 0x08

#define  TN93_MAX_DIST 1000.0

using namespace std;

StringBuffer names,
sequences;

Vector       nameLengths,
seqLengths,
workingNodes,
nodeParents;

VectorDouble distanceEstimates;

const  char ValidChars[]       = "ACGTURYSWKMBDHVN?-",
            ValidCharsAA[]     = "ACDEFGHIKLMNPQRSTVWYBZX?-";
            
unsigned char   * resolveTheseAmbigs = (unsigned char   *)calloc (256,sizeof (unsigned char));
                  
static char empty        [] = "";


const long   resolutions [][4] = { {1,0,0,0},
  {0,1,0,0},
  {0,0,1,0},
  {0,0,0,1},
  {0,0,0,1}, // U - 4
  {1,0,1,0}, //RESOLVE_A | RESOLVE_G, // R - 5
  {0,1,0,1}, //RESOLVE_C | RESOLVE_T, // Y - 6
  {0,1,1,0}, //RESOLVE_C | RESOLVE_G, // S - 7
  {1,0,0,1}, //RESOLVE_A | RESOLVE_T, // W - 8
  {0,0,1,1}, //RESOLVE_G | RESOLVE_T, // K - 9
  {1,1,0,0}, //RESOLVE_A | RESOLVE_C, // M - 10
  {0,1,1,1}, // RESOLVE_C | RESOLVE_G | RESOLVE_T, // B - 11
  {1,0,1,1}, //RESOLVE_A | RESOLVE_G | RESOLVE_T, // D - 12
  {1,1,0,1}, //RESOLVE_A | RESOLVE_C | RESOLVE_T, // H - 13
  {1,1,1,0}, // RESOLVE_A | RESOLVE_C | RESOLVE_G, // V - 14
  {1,1,1,1}, // RESOLVE_A | RESOLVE_C | RESOLVE_G | RESOLVE_T , // N - 15
  {1,1,1,1}, //RESOLVE_A | RESOLVE_C | RESOLVE_G | RESOLVE_T , // ? - 16
  {0,0,0,0} // GAP
};




/*A.................Ala.................Alanine
B.................Asx.................Aspartic acid or Asparagine
C.................Cys.................Cysteine
D.................Asp.................Aspartic Acid
E.................Glu.................Glutamic Acid
F.................Phe.................Phenylalanine
G.................Gly.................Glycine
H.................His.................Histidine
I.................Ile.................Isoleucine
K.................Lys.................Lysine
L.................Leu.................Leucine
M.................Met.................Methionine
N.................Asn.................Asparagine
P.................Pro.................Proline
Q.................Gln.................Glutamine
R.................Arg.................Arginine
S.................Ser.................Serine
T.................Thr.................Threonine
V.................Val.................Valine
W.................Trp.................Tryptophan
X.................Xaa.................Any amino acid
Y.................Tyr.................Tyrosine
Z.................Glx.................Glutamine or Glutamic acid*/
 
 
// ACDEFGHIKLMNPQRSTVWYBZ?- 

const long   resolutions_AA [][20] = { {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
//                                       A C D E F G H I L K M N P Q R S T V W Y
                                        {0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, // B
                                        {0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0}, // Z
                                        {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}, // X
                                        {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}, // ?
                                        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}  // -

};

#define N_CHAR 15
#define GAP    17
#define GAP_AA 24

const  double   resolutionsCount[] = { 1.f,
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
  0.f
};

const  double   resolutionsCount_AA[] = { 
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  1.f,
  0.5f,
  0.5f,
  0.05f,
  0.05f,
  0.f
};


char validFlags [256];

  //---------------------------------------------------------------

void initAlphabets (bool doAminoAcid, char * resolutionSubset) {
  for (int i = 0; i < 256; i++)
    validFlags [i] = -1;
  
  if (doAminoAcid) {
    for (unsigned int i = 0; i < strlen (ValidCharsAA); i++)
      validFlags [(unsigned char)ValidCharsAA[i]] = i;
    
    
  } else {  
     for (unsigned int i = 0; i < strlen (ValidChars); i++)
      validFlags [(unsigned char)ValidChars[i]] = i;
      
    if (resolutionSubset) {
      unsigned long subset_length = strlen (resolutionSubset);
      for (unsigned long rc = 0; rc < subset_length; rc++) {
        unsigned char rcc = toupper( (resolutionSubset[rc]));
        if (validFlags[rcc] > 3) {
          resolveTheseAmbigs[validFlags[rcc]] = 1;
        }
      }
    }
  }
  
}

  //---------------------------------------------------------------

const char unmap_char (unsigned char c, bool do_aa) {
  return do_aa ? ValidCharsAA[c] : ValidChars[c];
}

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
  static unsigned long mag01[2]= {0x0UL, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
  if (mti >= N) { /* generate N words at one time */
    int kk;
    
    if (mti == N+1)   /* if init_genrand() has not been called, */
      init_genrand(5489UL); /* a default initial seed is used */
    
    for (kk=0; kk<N-M; kk++) {
      y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
      mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk<N-1; kk++) {
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

const long* resolve_char (unsigned char c, bool do_aa, bool ignore_ambigs) {
  
  if (do_aa) {
    if (ignore_ambigs && c >= 20) {
      return resolutions_AA[GAP_AA];
    }
    return resolutions_AA [c];
  } 
  
  if (ignore_ambigs && c >= 4) {
    return resolutions[GAP];
  }
  
  return resolutions [c];
  
}

  //---------------------------------------------------------------

const double resolution_count (unsigned char c, bool do_aa) {
  
  if (do_aa) {
    return resolutionsCount_AA[c];
  } 
   return resolutionsCount [c];
  
}

  //---------------------------------------------------------------

long stringLength (Vector& lengths, unsigned long index)
{
  if (index < lengths.length() - 1)
    return lengths.value(index+1) - lengths.value(index) - 1;
  
  return -1;
}

  //---------------------------------------------------------------

char* stringText (const StringBuffer& strings, const Vector& lengths, unsigned long index)
{
  if (index < lengths.length() - 1)
    return strings.getString() + lengths.value(index);
  
  return empty;
}

  //---------------------------------------------------------------

void merge_two_sequences (const char* source, char* target, const long sequence_length) {
  for (long char_index = 0; char_index < sequence_length; char_index ++) {
    if (target[char_index] == GAP &&  source[char_index]!=GAP) {
      target[char_index] = source[char_index];
    }
  }
}


  //---------------------------------------------------------------

long perfect_match (const char* source, char* target, const long sequence_length) {
  long matched_bases = 0;
  for (long c = 0; c < sequence_length; c++) {
    char c1 = source[c],
    c2 = target[c];
    
    if (c1 == GAP || c2 == GAP) continue;
    
    if (c1 != c2) {
      return -1;
    }
    
    matched_bases ++;
  }
  return matched_bases;
}

/*---------------------------------------------------------------------------------------------------- */

double		computeTN93 (char * s1, char *s2,  unsigned long L, char matchMode, long* randomize, long min_overlap,
                       unsigned long* histogram, double slice, unsigned long hist_size, long count1, long count2) {
  
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
  
  for (unsigned long i = 0; i < 4UL; i++)
    for (unsigned long j = 0; j < 4UL; j++)
      pairwiseCounts[i][j] = 0.;
  
  
  for (unsigned long p = 0; p < L; p++) {
    unsigned char c1, c2;
    
    if (randomize == NULL) {
      c1 = s1[p];
      c2 = s2[p];
    } else {
      long pi = randomize[p];
      c1 = s1[pi];
      c2 = s2[pi];
    }
    
    if (c1 == GAP || c2 == GAP) {
      if (matchMode != GAPMM) {
        continue;
      } else {
        if (c1 == GAP && c2 == GAP)
          continue;
        else {
          if (c1 == GAP) {
            c1 = N_CHAR;
          } else {
            c2 = N_CHAR;
          }
        } 
      } 
    }
    
    if (c1 < 4) {
      if (c2 < 4) {
        pairwiseCounts [c1][c2] += 1.;
      }
      else {
        if (matchMode != SKIP) {
          if (resolutionsCount[c2] > 0.) {
            if (matchMode == RESOLVE || matchMode == SUBSET && resolveTheseAmbigs[c2])
              if (resolutions[c2][c1]) {
                pairwiseCounts[c1][c1] += 1.;
                continue;
              }
            
            if (resolutions[c2][0])
              pairwiseCounts[c1][0] += resolutionsCount[c2];
            if (resolutions[c2][1])
              pairwiseCounts[c1][1] += resolutionsCount[c2];
            if (resolutions[c2][2])
              pairwiseCounts[c1][2] += resolutionsCount[c2];
            if (resolutions[c2][3])
              pairwiseCounts[c1][3] += resolutionsCount[c2];
          }
        }
      }
    }
    else {
      if (matchMode != SKIP) {
        if (c2 < 4) {
          if (resolutionsCount[c1] > 0.) {
            if (matchMode == RESOLVE || matchMode == SUBSET && resolveTheseAmbigs[c1]) {
              if (resolutions[c1][c2]) {
                pairwiseCounts[c2][c2] += 1.;
                continue;
              }
            }
            
            if (resolutions[c1][0])
              pairwiseCounts[0][c2] += resolutionsCount[c1];
            if (resolutions[c1][1])
              pairwiseCounts[1][c2] += resolutionsCount[c1];
            if (resolutions[c1][2])
              pairwiseCounts[2][c2] += resolutionsCount[c1];
            if (resolutions[c1][3])
              pairwiseCounts[3][c2] += resolutionsCount[c1];
          }
       } else {
            // ambig and ambig
          double norm = resolutionsCount[c1] * resolutionsCount[c2];
            //cout << int(c1) << ":" << int(c2) << "/" << norm << endl;
          if (norm > 0.0) {
            if (matchMode == RESOLVE || matchMode == SUBSET && resolveTheseAmbigs[c1] && resolveTheseAmbigs[c2]) {
              long matched_count = 0L,
                   positive_match [4] = {0,0,0,0};
              for (long i = 0; i < 4L; i ++) {
                if (resolutions[c1][i] && resolutions[c2][i]) {
                  matched_count ++;
                  positive_match[i] = 1;
                }
              }
              
              if (matched_count > 0L) {
                double norm2 = 1./matched_count;
                
                for (long i = 0; i < 4L; i ++) {
                  if (positive_match[i]) {
                    pairwiseCounts[i][i] += norm2;
                  }
                }
                continue;
              }
            }
            
            for (long i = 0; i < 4L; i ++) {
              if (resolutions[c1][i]) {
                for (long j = 0; j < 4L; j ++) {
                  if (resolutions [c2][j]) {
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
  
  
  for (unsigned long c1 = 0; c1 < 4; c1++) {
    for (unsigned long c2 = 0; c2 < 4; c2++) {
      double pc = pairwiseCounts[c1][c2];
      totalNonGap   += pc;
      nucFreq [c1]  += pc;
      nucFreq [c2]  += pc;
    }
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
  
  
  double dist;
  
  
  if (useK2P)
      {
    ti = AG+CT;
    AG = 1.-2.*ti-tv;
    CT = 1.-2.*tv;
    if (AG > 0. && CT > 0.)
      dist = -0.5*log(AG)-0.25*log(CT);
    else  
      dist = TN93_MAX_DIST;
      }
  else
      {
    AG	= 1.-AG/K1 - 0.5*tv/fR;
    CT	= 1.-CT/K2 - 0.5*tv/fY;
    tv  = 1.-0.5 * tv/fY/fR;
    if (AG > 0. && CT > 0. && tv > 0)
      dist = - K1*log(AG) - K2*log(CT) - K3 * log (tv);
    else
      dist = TN93_MAX_DIST;
      }
  
  if (histogram) {
    long index = floor(dist * slice);
    if (index >= hist_size) {
      index = hist_size - 1;
    }
    histogram[index] += count1 * count2;
  }
  
  return dist;
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

int readFASTA (FILE* F, char& automatonState,  StringBuffer &names, 
               StringBuffer& sequences, Vector &nameLengths, Vector &seqLengths, 
               long& firstSequenceLength, bool oneByOne, 
               Vector* sequenceInstances, char sep,
               double include_prob) {
  
  unsigned long up_to = 0L;
  
  if (oneByOne) {
    sequences.resetString();
    names.resetString();
    include_prob = 1.;
  }
  
  if (include_prob < 1.) {
    up_to = RAND_RANGE * include_prob;
  }
  
  bool include_me = true;
  
  while (1) {
    int currentC = fgetc (F);
      //cout << "State: " << int(automatonState) << "/'" << char(currentC) << "'" << endl;
    if (feof (F))
      break;
    switch (automatonState) {
      case 0: {
        if (currentC == '>' || currentC == '#') {
          automatonState = 1;
          if (include_prob < 1.) {
            include_me = genrand_int32() < up_to;
          }
        }
        break;
      }
      case 1: {
        if (currentC == '\n' || currentC == '\r') {
          if (include_me) {
            names.appendChar   ('\0');
            long this_name_l;
            if (oneByOne) {
              this_name_l = names.length ()-1;        
              
            } else {
              nameLengths.appendValue (names.length());
              this_name_l = stringLength (nameLengths, nameLengths.length()-2);
            }
            
            if (this_name_l <= 0) {
              cerr << "Sequence names must be non-empty." << endl;
              return 1;
            }
            
            if (sequenceInstances) {
              unsigned long count = 1L;
              if (this_name_l >= 3) {
                long sep_loc = 0, ll = names.length();
                for (sep_loc = 2; sep_loc < this_name_l; sep_loc ++) {
                  if (names.getChar(ll-sep_loc-1) == sep) {
                    break;
                  }
                }
                if (sep_loc < this_name_l) {
                  count = atoi (names.getString() + (ll-sep_loc));
                }
                if (count < 1) {
                  count = 1;
                }
              }
              if (oneByOne) {
                sequenceInstances->resetVector();
              }
              sequenceInstances->appendValue(count);
            }
          }
          
          automatonState = 2;
        }
        else {
          if (include_me) {
            names.appendChar(currentC);
          }
        }
        break;
      }
      case 2: {
        currentC = toupper (currentC);
        if (validFlags [currentC] >= 0) {
          if (include_me)
            //cout << "Append " << currentC << endl;
            sequences.appendChar (validFlags [currentC]);
        }
        else {
          if (currentC == '>' || currentC == '#') {
            automatonState = 1;
            if (include_me) {
              if (oneByOne) {
                if (firstSequenceLength == 0) {
                  firstSequenceLength = sequences.length()-1;
                } 
                  //cerr << endl << "Returning a sequence" << endl;
                automatonState = 0;
                sequences.appendChar ('\0');
                ungetc (currentC, F);
                return 2;
              }
              addASequenceToList (sequences, seqLengths, firstSequenceLength, names, nameLengths);
            }
            if (include_prob < 1.) {
              include_me = genrand_int32() < up_to;
            }
          }
        }
        break;
      }
    }
    
      }
  
  if (automatonState == 2 || (oneByOne && automatonState == 0)) {
    if (!oneByOne && include_me) {
      addASequenceToList (sequences, seqLengths, firstSequenceLength, names, nameLengths);
    }
    automatonState = 1;
  } else {
    cerr << "Unexpected end of file: state " << int(automatonState) << endl;
    return 1;
  }
  return 0;
}

  //---------------------------------------------------------------

void dump_fasta (const char* mapped_characters, const long firstSequenceLength, FILE *output, bool newln, bool is_prot, unsigned long from, unsigned long to) {
  if (from > 0 || to > 0) {
    for (long c = from; c <= to; c++) {
      fputc( is_prot ? ValidCharsAA[mapped_characters[c]] : ValidChars[mapped_characters[c]], output);    
    }
  } else {
    for (long c = 0; c < firstSequenceLength; c++) {
      fputc( is_prot ? ValidCharsAA[mapped_characters[c]] : ValidChars[mapped_characters[c]], output);
    }
  }
  if (newln) {
    fprintf (output, "\n"); 
  }
}

  //---------------------------------------------------------------

void dump_sequence_fasta (unsigned long index, FILE* output, long firstSequenceLength, double * d, bool is_prot, unsigned long from, unsigned long to) {
  if (d) {
    fprintf (output, ">%s [%g, %g]\n", stringText (names, nameLengths, index), d[0], d[1]);
  } else {
    fprintf (output, ">%s\n", stringText (names, nameLengths, index));
  }
  
  dump_fasta (stringText (sequences, seqLengths, index), firstSequenceLength, output, true, is_prot, from, to);
}



