#ifndef SUBQUADRATIC_H
#define SUBQUADRATIC_H

#include "tn93_shared.h"
#include "argparse.hpp"

int run_subquadratic_tn93 (argparse::args_t &args,
                           StringBuffer &sequences,
                           Vector &seqLengths,
                           StringBuffer &names,
                           Vector &nameLengths,
                           sequence_gap_structure *sequence_descriptors,
                           long firstSequenceLength,
                           Vector &counts,
                           int resolutionOption,
                           unsigned long seqLengthInFile1,
                           unsigned long seqLengthInFile2);

#endif // SUBQUADRATIC_H
