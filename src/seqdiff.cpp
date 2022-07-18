
#include "tn93_shared.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#include "argparse_seqdiff.hpp"
#include <map>
#include <vector>
#include<algorithm>

using namespace std;
using namespace argparse;


struct cmpMatchPatterns {
    bool operator()(const Vector* a, const Vector* b) const {
        const long al = a->length ();
        if (al < b->length()) return true;
        if (al > b->length()) return false;
        for (long i = 0; i < al; i++) {
            long diff = a->value(i) - b->value (i);
            if (diff < 0) return true;
            if (diff > 0) return false;
       }
       return false;
    }
};

//---------------------------------------------------------------

int main(int argc, const char *argv[]) {
  args_t args = args_t(argc, argv);

  long firstSequenceLength = 0;

  char automatonState = 0;
  // 0 - between sequences
  // 1 - reading sequence name
  // 2 - reading sequence

  nameLengths.appendValue(0);
  seqLengths.appendValue(0);

  initAlphabets(false, NULL);
  if (readFASTA(args.reference, automatonState, names, sequences, nameLengths,
                  seqLengths, firstSequenceLength) == 1)
    return 1;
    
  automatonState = 0;
    
    
  if (readFASTA(args.input, automatonState, names, sequences, nameLengths,
                seqLengths, firstSequenceLength, false, NULL, ':', 1.0, !args.quiet) == 1)
    return 1;

  unsigned long sequenceCount = seqLengths.length() - 1;

  if (!args.quiet) {
      cerr << "Read " << sequenceCount - 1 << " sequences of length " << firstSequenceLength <<"nt from input file " << endl;
  }

  sequence_gap_structure *sequence_descriptors = NULL;

  double percentDone = 0.;

  sequence_descriptors = new sequence_gap_structure[sequenceCount];
  for (unsigned long sid = 0UL; sid < sequenceCount; sid++) {
    sequence_descriptors[sid] = describe_sequence(
        stringText(sequences, seqLengths, sid), firstSequenceLength);
  }

  if (!args.quiet)
    cerr << endl << "Progress: ";

  int resolutionOption;
  switch (args.ambig) {
  case resolve:
    resolutionOption = RESOLVE;
    break;
  case all:
    resolutionOption = MISMATCH;
    break;
  case informative:
    resolutionOption = INFORMATIVE;
    break;
  }

  time_t before, after;
  time(&before);

  const char *ref_name = stringText(names, nameLengths, 0),
             *ref_seq = stringText(sequences, seqLengths, 0);
    
  const sequence_gap_structure *ref_gaps = &sequence_descriptors[0];
    
  std::map <Vector*, long, cmpMatchPatterns> pattern_map;
  std::vector <long> assignments;
  std::vector <Vector*> patterns;
  std::vector <std::vector<long>> copies;
    
    
  for (long seq_id = 1L; seq_id < sequenceCount; seq_id ++) {
      
      if (!args.quiet && (seq_id * 100. / sequenceCount - percentDone > 0.1 ||
                          seq_id == (long)sequenceCount - 1)) {
        {
          time(&after);
          percentDone = seq_id * 100. / sequenceCount;
          cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress"
                  ":"
               << setw(8) << percentDone << "% (" << setw(12) << std::setprecision(3)
               << seq_id / difftime(after, before) << " seqs/sec)";

          after = before;
        }
      }
      
      Vector* differences = new Vector;
      const char *seq_data = stringText(sequences, seqLengths, seq_id);

      computeDifferences(ref_seq, seq_data, firstSequenceLength, resolutionOption, *differences, ref_gaps, &sequence_descriptors[seq_id]);
      
      auto intert_result = pattern_map.insert({differences, pattern_map.size()});
      bool did_insert = std::get<1>(intert_result);
      long my_code;
      if (!did_insert) {
          auto value = std::get<0>(intert_result);
          my_code = std::get<1>(*value);
          assignments.push_back(my_code);
          delete differences;
      } else {
          my_code = pattern_map.size()-1;
          assignments.push_back(pattern_map.size()-1);
          patterns.push_back (differences);
      }
      if (my_code >= copies.size()) {
          std::vector<long> new_vector;
          copies.push_back (new_vector);
      }
      copies[my_code].push_back(seq_id);
  }
    
  
    
  if (args.out_format == json) {
      if (args.compress) {
          std::map <std::string, long> differences_map;
          std::vector <std::string> unique_differences;
          std::vector <std::vector<long>> recoded_patterns;
          
          for (Vector* p : patterns) {
              std::vector <long> recoded;
              for (long i = 0; i < p->length(); i++) {
                  std::string mut;
                  long     position;
                  unsigned character;
                  unpack_difference(p->value (i), position, character);
                  mut += to_string(position+1);
                  mut += unmap_char(character);
                  auto intert_result = differences_map.insert({mut, differences_map.size()});
                  bool did_insert = std::get<1>(intert_result);
                  long index;
                  if (did_insert) {
                      unique_differences.push_back (mut);
                      //cout << mut << endl;
                      index = differences_map.size() - 1;
                  } else {
                      auto value = std::get<0>(intert_result);
                      index = std::get<1>(*value);
                  }
                  recoded.push_back (index);
              }
              recoded_patterns.push_back (recoded);
          }
          
          // write out mutational codes
          fprintf (args.output, "{\"reference\":\"");
          dump_fasta (ref_seq, firstSequenceLength, args.output, false, false,0L,0L);
          fprintf (args.output, "\",\"mutations\" : [");
          auto mutation_string = unique_differences.begin();
          fprintf (args.output, "\"%s\"", (*mutation_string).c_str());
          mutation_string ++;
          std::for_each(mutation_string, unique_differences.end(), [&](const std::string& s) {fprintf (args.output,",\"%s\"",s.c_str()); });
          fprintf (args.output, "], \"patterns\" : [");
          bool comma = false;
          std::for_each(recoded_patterns.begin(), recoded_patterns.end(),[&](const std::vector<long>& pattern) {
              if (comma) fprintf (args.output, ",");
              long comma_inner = false;
              fprintf (args.output, "[");
              std::for_each(pattern.begin(), pattern.end(), [&](const long mi) {
                  if (comma_inner) fprintf (args.output, ",");
                  fprintf (args.output, "%ld", mi);
                  comma_inner = true;
              });
              fprintf (args.output, "]");
              comma = true;
          });
          fprintf (args.output, "], \"sequences\" : [");
          comma = false;
          std::for_each(copies.begin(), copies.end(),[&](const std::vector<long>& pattern) {
              if (comma) fprintf (args.output, ",");
              long comma_inner = false;
              fprintf (args.output, "[");
              std::for_each(pattern.begin(), pattern.end(), [&](const long mi) {
                  if (comma_inner) fprintf (args.output, ",");
                  fprintf (args.output, "\"%s\"", stringText(names, nameLengths, mi));
                  comma_inner = true;
              });
              fprintf (args.output, "]");
              comma = true;
          });
          fprintf (args.output, "]}");
      }
      
      /*for (std::vector<long> c: copies ) {
          cout << c.size() << endl;
      }*/
      
      
  }
    
  delete [] sequence_descriptors;
  return 0;
}
