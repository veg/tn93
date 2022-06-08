#ifndef 	__TN93SHARED__
#define 	__TN93SHARED__

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <climits>
#include "stringBuffer.h"

using namespace std;

#define  RESOLVE_A      0x01
#define  RESOLVE_C      0x02
#define  RESOLVE_G      0x04
#define  RESOLVE_T      0x08


#define  RESOLVE        0
#define  AVERAGE        1
#define  SKIP           2
#define  GAPMM          3
#define  SUBSET         4
#define  MISMATCH       5
#define  INFORMATIVE    6

#define RAND_RANGE 0xffffffffUL /* Maximum value returned by genrand_int32 */

#define MIN(a,b) (a) < (b) ? (a) : (b)
#define MAX(a,b) (a) > (b) ? (a) : (b)

struct sequence_gap_structure {
  
  long first_nongap,
       last_nongap,
       resolved_start,
       resolved_end;
  
  sequence_gap_structure (void) {
    first_nongap    = LONG_MAX;
    last_nongap     = 0L;
    resolved_start = 0L;
    resolved_end   = 0L;
  }
  
};

void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
double		computeTN93 (const char * s1, const char *s2,  const unsigned long L, const char matchMode, const long * randomize, const long min_overlap, unsigned long* = NULL, const double = 0.0, const unsigned long cnt = 0, const long count1 = 1, const long count2 = 1, const sequence_gap_structure * = NULL, const sequence_gap_structure * = NULL);

long   computeDifferences (const char * s1,
                           const char *s2,
                           const unsigned long L,
                           const char matchMode,
                           Vector& storage,
                           const sequence_gap_structure * = NULL,
                           const sequence_gap_structure * = NULL);


long stringLength (Vector& lengths, unsigned long index);
char* stringText (const StringBuffer& strings, const Vector& lengths, unsigned long index);
void addASequenceToList (StringBuffer& sequences, Vector& seqLengths, long &firstSequenceLength, StringBuffer& names, Vector& nameLengths);
int readFASTA (FILE* F, char& automatonState,  StringBuffer &names, StringBuffer& sequences, Vector &nameLengths, Vector &seqLengths, long& firstSequenceLength, bool oneByOne = false,  Vector* sequenceInstances = NULL, char sep = ':', double include_prob = 1.0, bool show_progress = false);
void dump_sequence_fasta (unsigned long index, FILE* output, long firstSequenceLength, double * d = NULL, bool = false, unsigned long from = 0L, unsigned long to = 0L);
void initAlphabets(bool = false, char * = NULL, bool id_map = false);
void merge_two_sequences (const char* source, char* target, const long sequence_length);
long perfect_match (const char* source, char* target, const long sequence_length);
void dump_fasta (const char* source, const long sequence_length, FILE* output, bool newln = true, bool = false, unsigned long from = 0L, unsigned long to = 0L);

struct sequence_gap_structure describe_sequence (const char* source, const unsigned long sequence_length, const unsigned long char_count = 4UL);

const long * resolve_char (unsigned char, bool = false, bool = true);
const double resolution_count (unsigned char, bool = false);
const char unmap_char (unsigned char, bool = false);
inline void unpack_difference (long diff, long& location, unsigned& alt) {
    location = diff >> 8;
    alt = diff & 0xff;
}


extern StringBuffer names,
       sequences;

extern unsigned char * resolveTheseAmbigs;

extern double   resolve_fraction;

extern Vector       nameLengths,
       seqLengths,
       workingNodes,
       nodeParents;

extern VectorDouble distanceEstimates;
extern const  double  resolutionsCount [];


#endif
