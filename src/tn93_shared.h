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

#define  RESOLVE_A      0x01
#define  RESOLVE_C      0x02
#define  RESOLVE_G      0x04
#define  RESOLVE_T      0x08


#define  RESOLVE        0
#define  AVERAGE        1
#define  SKIP           2
#define  GAPMM          3

#define RAND_RANGE 0xffffffffUL /* Maximum value returned by genrand_int32 */

void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
double		computeTN93 (char * s1, char *s2,  unsigned long L, char matchMode, long* randomize, long min_overlap);
long stringLength (Vector& lengths, unsigned long index);
char* stringText (StringBuffer& strings, Vector& lengths, unsigned long index);
void addASequenceToList (StringBuffer& sequences, Vector& seqLengths, long &firstSequenceLength, StringBuffer& names, Vector& nameLengths);
int readFASTA (FILE* F, char& automatonState,  StringBuffer &names, StringBuffer& sequences, Vector &nameLengths, Vector &seqLengths, long& firstSequenceLength, bool oneByOne = false);
void dump_sequence_fasta (unsigned long index, FILE* output, long firstSequenceLength, double * d = NULL);
void initAlphabets(void);
void merge_two_sequences (const char* source, char* target, const long sequence_length);
long perfect_match (const char* source, char* target, const long sequence_length);
void dump_fasta (const char* source, const long sequence_length, FILE* output, bool newln = true);

extern StringBuffer names,
       sequences;

extern Vector       nameLengths,
       seqLengths,
       workingNodes,
       nodeParents;

extern VectorDouble distanceEstimates;

#endif