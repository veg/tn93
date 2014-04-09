#ifndef 	__TN93SHARED__
#define 	__TN93SHARED__

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <unistd.h>
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

#define RAND_RANGE 0xffffffffUL /* Maximum value returned by genrand_int32 */

void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
double		computeTN93 (char * s1, char *s2,  unsigned long L, char matchMode, long* randomize, long min_overlap, unsigned long* = NULL, double = 0.0, unsigned long cnt = 0, long count1 = 1, long count2 = 1);
long stringLength (Vector& lengths, unsigned long index);
char* stringText (const StringBuffer& strings, const Vector& lengths, unsigned long index);
void addASequenceToList (StringBuffer& sequences, Vector& seqLengths, long &firstSequenceLength, StringBuffer& names, Vector& nameLengths);
int readFASTA (FILE* F, char& automatonState,  StringBuffer &names, StringBuffer& sequences, Vector &nameLengths, Vector &seqLengths, long& firstSequenceLength, bool oneByOne = false,  Vector* sequenceInstances = NULL, char sep = ':', double include_prob = 1.0);
void dump_sequence_fasta (unsigned long index, FILE* output, long firstSequenceLength, double * d = NULL, bool = false, unsigned long from = 0L, unsigned long to = 0L);
void initAlphabets(bool = false, char * = NULL);
void merge_two_sequences (const char* source, char* target, const long sequence_length);
long perfect_match (const char* source, char* target, const long sequence_length);
void dump_fasta (const char* source, const long sequence_length, FILE* output, bool newln = true, bool = false, unsigned long from = 0L, unsigned long to = 0L);

const long * resolve_char (unsigned char, bool = false, bool = true);
const double resolution_count (unsigned char, bool = false);
const char unmap_char (unsigned char, bool = false);

extern StringBuffer names,
       sequences;

extern unsigned char * resolveTheseAmbigs;

extern Vector       nameLengths,
       seqLengths,
       workingNodes,
       nodeParents;

extern VectorDouble distanceEstimates;


#endif
