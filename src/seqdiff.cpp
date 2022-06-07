
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
                seqLengths, firstSequenceLength) == 1)
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
               << setw(8) << percentDone << "% " << setw(12) << std::setprecision(3)
               << seq_id / difftime(after, before) << " evals/sec)";

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
          fprintf (args.output, "{\"mutations\" : [");
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
  //fclose (args.input);
  //fclose (args.output);
  return 0;

  /*if (args.cluster_type == all) {

    auto outer_iterator = remaining.begin();

    while (outer_iterator != remaining.end()) {

      unsigned long seq1 = *outer_iterator;

      char *base_sequence = stringText(sequences, seqLengths, seq1);

      Vector *current_cluster = clusters.append_a_cluster(seq1);

      outer_iterator = remaining.erase(outer_iterator);

      auto inner_iterator = remaining.begin();

      while (inner_iterator != remaining.end()) {

        unsigned long seq2 = *inner_iterator;
        char *test_sequence = stringText(sequences, seqLengths, seq2);

        try {
          double distance = computeTN93(
              base_sequence, test_sequence, firstSequenceLength,
              resolutionOption, NULL, args.overlap, NULL, HISTOGRAM_SLICE,
              HISTOGRAM_BINS, 1L, 1L, &sequence_descriptors[seq1],
              &sequence_descriptors[seq2]);

          if (distance <= args.distance) {
            // check if this can be merged into the current cluster
            for (unsigned long current_cluster_seq = 1UL;
                 current_cluster_seq < current_cluster->length();
                 current_cluster_seq++) {
              unsigned long csi = current_cluster->value(current_cluster_seq);
              double d = computeTN93(
                  test_sequence, stringText(sequences, seqLengths, csi),
                  firstSequenceLength, resolutionOption, NULL, args.overlap,
                  NULL, HISTOGRAM_SLICE, HISTOGRAM_BINS, 1L, 1L,
                  &sequence_descriptors[seq2], &sequence_descriptors[csi]);
                if (d > args.distance || d < 0.) {
                    throw(current_cluster_seq);
                }
            }
          } else {
            throw(0UL);
          }
          current_cluster->appendValue(seq2);
          //cerr << "removing " << seq2 << " : " << remaining.size() << " / " << current_cluster->length() << endl;
          inner_iterator = remaining.erase(inner_iterator);
            if (!args.quiet) {
              time(&after);
              double dt = difftime(after, before);
              if (dt > 0.) {
                percentDone =
                    (sequenceCount - remaining.size()) * 100. / sequenceCount;
                cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                        "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                        "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress"
                        ":"
                     << setw(8) << percentDone << "% (" << setw(8) << clusters.size()
                     << " clusters created, " << setw(12) << std::setprecision(3)
                     << ((double)sequenceCount - remaining.size()) / dt
                     << " sequences clustered/sec)";

                after = before;
              }
            }
        } catch (unsigned long e) {
          //cerr << "mismatch at " << e << endl;
          inner_iterator++;
        }
      }

      outer_iterator = remaining.begin();
      if (!args.quiet) {
        time(&after);
        double dt = difftime(after, before);
        if (dt > 0.) {
          percentDone =
              (sequenceCount - remaining.size()) * 100. / sequenceCount;
          cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress"
                  ":"
               << setw(8) << percentDone << "% (" << setw(8) << clusters.size()
               << " clusters created, " << setw(12) << std::setprecision(3)
               << ((double)sequenceCount - remaining.size()) / dt
               << " sequences clustered/sec)";

          after = before;
        }
      }
    }
  } else {
    auto outer_iterator = remaining.begin();

    // the idea here is to first iterate over all existing clusters
    // if the remaining sequence can be joined to some of the clusters,
    // then these clusters will be joined together and augmented with the new
    // sequence

    // failing that, the new sequence gets its own cluster

    unsigned long mergers = 0UL;

    while (outer_iterator != remaining.end()) {

      unsigned long seq1 = *outer_iterator;
      char *base_sequence = stringText(sequences, seqLengths, seq1);

      clusters.remove_blanks();

      set<unsigned long> join_to;

      outer_iterator = remaining.erase(outer_iterator);

#pragma omp parallel shared(clusters, sequence_descriptors, resolutionOption,  \
                            sequences, seqLengths, sequenceCount,              \
                            firstSequenceLength, args, nameLengths, names,     \
                            percentDone, cerr)
      for (unsigned long cc = 0UL; cc < clusters.size(); cc++) {
        Vector const *this_cluster = clusters[cc];

        try {
          for (unsigned long seq = 0UL; seq < this_cluster->length(); seq++) {
            unsigned long csi = this_cluster->value(seq);
            double d = computeTN93(
                base_sequence, stringText(sequences, seqLengths, csi),
                firstSequenceLength, resolutionOption, NULL, args.overlap, NULL,
                HISTOGRAM_SLICE, HISTOGRAM_BINS, 1L, 1L,
                &sequence_descriptors[seq1], &sequence_descriptors[csi]);
            if (d >= 0. && d <= args.distance) {
              throw(0);
            }
          }
        } catch (int e) {
#pragma omp critical
          join_to.insert(cc);
        }
      }

      if (join_to.empty()) {
        // create new cluster
        clusters.append_a_cluster(seq1);

      } else {
        auto cluster_iterator = join_to.begin();
        unsigned long first = *cluster_iterator;
        clusters[first]->appendValue(seq1);
        cluster_iterator++;
        while (cluster_iterator != join_to.end()) {
          // cout << "Merging " << first << " and " << *cluster_iterator << " on
          // " << seq1 << endl;
          clusters.merge_clusters(first, *cluster_iterator);
          mergers++;
          cluster_iterator++;
        }
      }

      outer_iterator = remaining.begin();
      if (!args.quiet) {
        time(&after);
        double dt = difftime(after, before);
        if (dt > 0.) {
          percentDone =
              (sequenceCount - remaining.size()) * 100. / sequenceCount;
          cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\bProgress:"
               << setw(8) << percentDone << "% (" << setw(8) << clusters.size()
               << "[" << setw(8) << mergers << "] clusters created [merged], "
               << setw(12) << std::setprecision(3)
               << ((double)sequenceCount - remaining.size()) / dt
               << " sequences clustered/sec)";

          after = before;
        }
      }
    }
  }

  if (!args.quiet) {
    cerr << endl << "Created " << clusters.size() << " clusters" << endl;
  }
  clusters.remove_blanks();

  if (args.output_mode == json) {
    FILE *print_to;

    if (args.trunk_path) {
      print_to = fopen(args.trunk_path, "w");
      if (!print_to) {
        cerr << "Failed to open '" << args.trunk_path << "' for writing"
             << endl;
        exit(0);
      }
    } else {
      print_to = stdout;
    }

    fprintf(print_to, "[");
    for (unsigned long cid = 0; cid < clusters.size(); cid++) {
      Vector *this_cluster = clusters[cid];
      if (cid) {
        fprintf(print_to, ",");
      }

      unsigned long max_nongap = 0UL, max_id = 0UL;

      fprintf(print_to, "\n\t{\n\t\t\"size\" : %ld, \n\t\t\"members\" : [",
              this_cluster->length());
      for (unsigned long p = 0; p < this_cluster->length(); p++) {
        fprintf(print_to, p ? ",\"%s\"" : "\"%s\"",
                stringText(names, nameLengths, this_cluster->value(p)));
        unsigned char *s = (unsigned char *)stringText(sequences, seqLengths,
                                                       this_cluster->value(p));
        unsigned long my_non_gap = 0UL;
        for (unsigned long ci = 0; ci < firstSequenceLength; ci++) {
          if (resolutionsCount[s[ci]] == 1.f) {
            my_non_gap++;
          }
        }
        if (my_non_gap > max_nongap) {
          max_nongap = my_non_gap;
          max_id = this_cluster->value(p);
        }
      }
      fprintf(print_to, "],\n\t\t\"centroid\" : \">%s\\n",
              stringText(names, nameLengths, max_id));
      dump_fasta(stringText(sequences, seqLengths, max_id), firstSequenceLength,
                 print_to, false);
      // dump_sequence_fasta (max_id, print_to, firstSequenceLength);
      fprintf(print_to, "\\n\"\n\t}");
    }
    fprintf(print_to, "\n]\n");

    if (args.trunk_path) {
      fclose(print_to);
    }

  } else {
    char *path_buffer = NULL;
    unsigned long pl = 0UL;
    if (args.trunk_path) {
      pl = strlen(args.trunk_path) + 128;
      path_buffer = new char[pl];
    }
    for (unsigned long cid = 0; cid < clusters.size(); cid++) {
      Vector *this_cluster = clusters[cid];

      FILE *print_to;

      if (args.trunk_path) {
        snprintf(path_buffer, pl, "%s.%lu", args.trunk_path, cid);
        // cout << path_buffer << endl;
        print_to = fopen(path_buffer, "w");
        if (!print_to) {
          cerr << "Failed to open '" << path_buffer << "' for writing" << endl;
          exit(0);
        }
      } else {
        print_to = stdout;
        cout << endl << "============================" << endl;
      }
      for (unsigned long p = 0; p < this_cluster->length(); p++) {
        long mapped_id = this_cluster->value(p);
        // cout << mapped_id << endl;
        dump_sequence_fasta(mapped_id, print_to, firstSequenceLength);
      }

      if (args.trunk_path) {
        fclose(print_to);
      }
    }
    if (args.trunk_path) {
      delete[] path_buffer;
    }
  }*/

  delete[] sequence_descriptors;
  return 0;
}
