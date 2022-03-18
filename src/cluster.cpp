
#include "tn93_shared.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#include "argparse_cluster.hpp"
#include <set>

using namespace std;
using namespace argparse;

#define HISTOGRAM_BINS 200
#define HISTOGRAM_SLICE ((double)HISTOGRAM_BINS)

template <typename datatype> void _swap(datatype &a, datatype &b) {
  datatype t = b;
  b = a;
  a = t;
}

//---------------------------------------------------------------

class Clusters {
public:
  Clusters(void) {
    allocated_clusters = 128UL;
    cluster_count = 0UL;
    unused = 0UL;
    list_of_clusters = new Vector *[allocated_clusters];
    for (unsigned long k = 0UL; k < allocated_clusters; k++) {
      list_of_clusters[k] = NULL;
    }
  }
  ~Clusters(void) {
    for (unsigned long k = 0UL; k < allocated_clusters; k++) {
      if (list_of_clusters[k]) {
        delete list_of_clusters[k];
      }
    }
    delete list_of_clusters;
  }

  unsigned long size(void) const { return cluster_count; }
  void remove_blanks(void) {
    if (unused > 0UL) {
      unsigned long shift = 0UL;
      for (unsigned long k = 0UL; k < allocated_clusters && shift <= unused;
           k++) {
        if (list_of_clusters[k] == NULL) {
          shift++;
        } else {
          _swap(list_of_clusters[k - shift], list_of_clusters[k]);
        }
      }
      unused = 0UL;
    }
  }

  void remove_a_cluster(unsigned long cluster_id) {
    if (cluster_id < allocated_clusters && list_of_clusters[cluster_id]) {
      unused++;
      delete list_of_clusters[cluster_id];
      list_of_clusters[cluster_id] = NULL;
    } else {
      cerr << "Requested the removal of an invalid cluster " << cluster_id
           << " out of " << allocated_clusters << " allocated clusters "
           << endl;
      exit(1);
    }
  }

  void merge_clusters(unsigned long id1, unsigned long id2) {
    if (id1 != id2 && id1 < allocated_clusters && id2 < allocated_clusters &&
        list_of_clusters[id1] && list_of_clusters[id2]) {
      list_of_clusters[id1]->appendVector(*list_of_clusters[id2]);
      list_of_clusters[id1]->sort();
      delete list_of_clusters[id2];
      list_of_clusters[id2] = NULL;
      cluster_count--;
      unused++;
    } else {
      cerr << "Requested the merger using an invalid cluster " << id1 << " and "
           << id2 << endl;
      exit(1);
    }
  }

  Vector *operator[](unsigned long index) {
    if (index < allocated_clusters) {
      return list_of_clusters[index];
    } else {
      cerr << "Requested an invalid cluster index " << index << " out of "
           << allocated_clusters << " allocated clusters " << endl;
      ;
      exit(1);
    }
  }

  Vector *append_a_cluster(unsigned long sequence_id) {
    /** allocate a new cluster and seed it with sequence sequence_id

     @param sequence_id

     return the new cluster vector
     */

    unsigned long slot = 0UL;

    if (unused > 0UL) {
      for (unsigned long k = 0UL; k < allocated_clusters; k++) {
        if (list_of_clusters[k] == NULL) {
          slot = k;
          break;
        }
        unused--;
      }
    } else {
      slot = cluster_count;
      if (cluster_count == allocated_clusters) {
        unsigned long new_allocation =
            allocated_clusters +
            (allocated_clusters > 512UL ? (allocated_clusters >> 2) : 128UL);
        Vector **new_vectors = new Vector *[new_allocation];

        for (unsigned long k = 0UL; k < allocated_clusters; k++) {
          new_vectors[k] = list_of_clusters[k];
        }
        for (unsigned long k = allocated_clusters; k < new_allocation; k++) {
          new_vectors[k] = NULL;
        }

        delete[] list_of_clusters;
        allocated_clusters = new_allocation;
        list_of_clusters = new_vectors;
      }
    }

    Vector *new_vector = (list_of_clusters[slot] = new Vector);
    new_vector->appendValue(sequence_id);
    cluster_count++;
    return new_vector;
  }

private:
  Vector **list_of_clusters;
  unsigned long allocated_clusters, cluster_count, unused;
};

//---------------------------------------------------------------

int compare_ul(const void *v1, const void *v2) {
  const unsigned long *ul1 = (const unsigned long *)v1,
                      *ul2 = (const unsigned long *)v2;

  if (*ul1 < *ul2) {
    return -1;
  }
  if (*ul1 > *ul2) {
    return 1;
  }
  return 0;
}

//---------------------------------------------------------------

int main(int argc, const char *argv[]) {
  args_t args = args_t(argc, argv);

  long firstSequenceLength = 0;

  set<unsigned long> remaining; // sequence ids not yet clustered

  Clusters clusters;

  char automatonState = 0;
  // 0 - between sequences
  // 1 - reading sequence name
  // 2 - reading sequence

  nameLengths.appendValue(0);
  seqLengths.appendValue(0);

  initAlphabets(false, args.ambigs_to_resolve);
  if (readFASTA(args.input, automatonState, names, sequences, nameLengths,
                seqLengths, firstSequenceLength) == 1)
    return 1;

  unsigned long sequenceCount = seqLengths.length() - 1,
                pairwise = (sequenceCount - 1L) * (sequenceCount) / 2;

  if (!args.quiet) {
    cerr << "Read " << sequenceCount << " sequences from input file " << endl
         << "Will perform approximately " << pairwise
         << " pairwise distance calculations";
  }

  sequence_gap_structure *sequence_descriptors = NULL;

  resolve_fraction = args.resolve_fraction;

  double percentDone = 0.;

  sequence_descriptors = new sequence_gap_structure[sequenceCount];
  for (unsigned long sid = 0UL; sid < sequenceCount; sid++) {
    sequence_descriptors[sid] = describe_sequence(
        stringText(sequences, seqLengths, sid), firstSequenceLength);
    if (args.first_regular || sid > 0UL) {
      remaining.insert(sid);
    }
  }

  if (!args.quiet)
    cerr << endl << "Progress: ";

  int resolutionOption;
  switch (args.ambig) {
  case resolve:
    resolutionOption = RESOLVE;
    break;
  case average:
    resolutionOption = AVERAGE;
    break;
  case skip:
    resolutionOption = SKIP;
    break;
  case gapmm:
    resolutionOption = GAPMM;
    break;
  case subset:
    resolutionOption = SUBSET;
    break;
  }

  time_t before, after;
  time(&before);

  if (args.cluster_type == all) {

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
  }

  delete[] sequence_descriptors;
  return 0;
}
