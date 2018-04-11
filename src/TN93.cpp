
#include "tn93_shared.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#include "argparse.hpp"

using namespace std;
using namespace argparse;

#define HISTOGRAM_BINS 200
#define HISTOGRAM_SLICE ((double)HISTOGRAM_BINS)

char sep = ':';

//---------------------------------------------------------------

void dump_histogram(ostream *outStream, const char *tag, unsigned long *hist) {
  if (tag) {
    (*outStream) << "\t\"Histogram " << tag << "\" : [";
  } else {
    (*outStream) << "\t\"Histogram\" : [";
  }
  for (unsigned long k = 0; k < HISTOGRAM_BINS; k++) {
    if (k) {
      (*outStream) << ',';
    }
    (*outStream) << '[' << (k + 1.) / HISTOGRAM_BINS << ',' << hist[k] << ']';
  }
  (*outStream) << "]" << endl;
}

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

  StringBuffer names, sequences;

  Vector nameLengths, seqLengths, counts, randomized_idx, recounts;

  bool randomized_fst = false;

  char automatonState = 0;
  // 0 - between sequences
  // 1 - reading sequence name
  // 2 - reading sequence

  init_genrand(time(NULL) + getpid());

  nameLengths.appendValue(0);
  seqLengths.appendValue(0);

  initAlphabets(false, args.ambigs_to_resolve);
  if (readFASTA(args.input1, automatonState, names, sequences, nameLengths,
                seqLengths, firstSequenceLength, false, &counts,
                args.counts_in_name, args.include_prob) == 1)
    return 1;

  unsigned long seqLengthInFile1 = seqLengths.length() - 1,
                seqLengthInFile2 = 0;

  if (args.input2) {
    automatonState = 0;
    if (readFASTA(args.input2, automatonState, names, sequences, nameLengths,
                  seqLengths, firstSequenceLength, false, &counts,
                  args.counts_in_name, args.include_prob) == 1)
      return 1;
    seqLengthInFile2 = seqLengths.length() - 1 - seqLengthInFile1;
  }

  bool do_fst = args.input2 && args.do_fst;
  bool report_self =
      args.input2 == NULL && args.report_self && args.do_count == false;

  unsigned long sequenceCount = seqLengths.length() - 1,
                pairwise = (args.input2 && !do_fst)
                               ? seqLengthInFile2 * seqLengthInFile1
                               : (sequenceCount - 1) * (sequenceCount) / 2;

  if (!args.quiet) {
    if (args.input2 == NULL) {
      cerr << "Read " << sequenceCount << " sequences of length "
           << firstSequenceLength << endl
           << "Will perform " << pairwise << " pairwise distance calculations";
    } else {
      cerr << "Read " << seqLengthInFile1 << " sequences from file 1 and "
           << seqLengthInFile2 << " sequences from file 2 of length "
           << firstSequenceLength << endl
           << "Will perform " << pairwise << " pairwise distance calculations";
      if (do_fst) {
        cerr << ". Running in F_ST mode";
      }
    }
  }

  sequence_gap_structure *sequence_descriptors = NULL;

  resolve_fraction = args.resolve_fraction;

  double percentDone = 0., max[3] = {0.0, 0.0, 0.0}, mean[3] = {0.0, 0.0, 0.0};

  long *randFlag = NULL;
  long *randSeqs = NULL;

  if (args.do_bootstrap) {
    if (args.input2 == NULL || args.do_bootstrap_two_files) {
      if (!args.quiet)
        cerr << endl << "Randomizing site order..." << endl;
      randFlag = new long[firstSequenceLength];
      bool *included = new bool[firstSequenceLength];
      for (long k = 0; k < firstSequenceLength; k++) {
        included[k] = false;
      }

      long normalizer = RAND_RANGE / firstSequenceLength, max_site = 0;

      for (long i = 0; i < firstSequenceLength; i++) {
        randFlag[i] = genrand_int32() / normalizer;
        if (!args.quiet)
          cerr << randFlag[i] << " ";
        included[randFlag[i]] = true;
        if (randFlag[i] > max_site)
          max_site = randFlag[i];
      }

      long total_resampled = 0;
      for (long k = 0; k < firstSequenceLength; k++) {
        total_resampled += included[k];
      }

      delete[] included;
      if (!args.quiet)
        cerr << endl
             << "Unique sites included in the resampled order "
             << total_resampled << endl;
    } else {
      randomized_fst = true;

      unsigned long total_count = 0UL, total_count_in_first = 0UL;

      for (unsigned long k = 0; k < sequenceCount; k++) {
        unsigned long my_count = counts.value(k);
        if (k < seqLengthInFile1) {
          total_count_in_first += my_count;
        }
        total_count += my_count;
      }

      unsigned long *expanded_counts = new unsigned long[total_count],
                    expanded_index = 0UL;

      for (unsigned long k = 0; k < sequenceCount; k++) {
        unsigned long my_count = counts.value(k);
        for (unsigned long i = 0; i < my_count; i++) {
          expanded_counts[expanded_index++] = k;
        }
      }

      for (unsigned long i = 0; i < total_count; i++) {
        long id = genrand_int32() % (total_count - i),
             t = expanded_counts[i + id];
        expanded_counts[i + id] = expanded_counts[i];
        expanded_counts[i] = t;
      }

      // void qsort (void* base, size_t num, size_t size,
      // int (*compar)(const void*,const void*))

      qsort(expanded_counts, total_count_in_first, sizeof(unsigned long),
            compare_ul);
      qsort(expanded_counts + total_count_in_first,
            total_count - total_count_in_first, sizeof(unsigned long),
            compare_ul);

      unsigned long last_sequence = expanded_counts[0], last_index = 0UL;

      for (unsigned long i = 0UL; i < total_count_in_first; i++) {
        if (expanded_counts[i] != last_sequence) {
          recounts.appendValue(i - last_index);
          last_index = i;
          randomized_idx.appendValue(last_sequence);
          last_sequence = expanded_counts[i];
        }
      }

      randomized_idx.appendValue(last_sequence);
      recounts.appendValue(total_count_in_first - last_index);
      seqLengthInFile1 = randomized_idx.length();

      last_sequence = expanded_counts[total_count_in_first];
      last_index = total_count_in_first;
      for (unsigned long i = total_count_in_first; i < total_count; i++) {
        if (expanded_counts[i] != last_sequence) {
          recounts.appendValue(i - last_index);
          last_index = i;
          randomized_idx.appendValue(last_sequence);
          // cerr << randomized_idx2.value (randomized_idx2.length ()-1) << "/"
          // << recount_2.value (recount_2.length()-1) << " (" << counts.value
          // (last_sequence) << ")" << endl;
          last_sequence = expanded_counts[i];
        }
      }

      randomized_idx.appendValue(last_sequence);
      recounts.appendValue(total_count - last_index);

      seqLengthInFile2 = randomized_idx.length() - seqLengthInFile1;

      if (!args.quiet)
        cerr << endl
             << "Randomized sequence assignment to datasets "
             << seqLengthInFile1 << " variants (" << total_count_in_first
             << ") to the first, and " << seqLengthInFile2 << " variants ("
             << (total_count - total_count_in_first) << ") to the second"
             << endl;

      sequence_descriptors = new sequence_gap_structure[sequenceCount];
      for (long sid = 0; sid < sequenceCount; sid++) {
        sequence_descriptors[sid] = describe_sequence(
            stringText(sequences, seqLengths, sid), firstSequenceLength);
      }

      sequenceCount = seqLengthInFile1 + seqLengthInFile2;
      pairwise = (sequenceCount - 1) * (sequenceCount) / 2;

      delete[] expanded_counts;
    }
  } else {
    sequence_descriptors = new sequence_gap_structure[sequenceCount];
    for (long sid = 0; sid < sequenceCount; sid++) {
      sequence_descriptors[sid] = describe_sequence(
          stringText(sequences, seqLengths, sid), firstSequenceLength);
    }
  }

  if (!args.quiet)
    cerr << endl << "Progress: ";

  long pairIndex = 0, foundLinks = 0;

  double *distanceMatrix = NULL;

  if (args.format != hyphy) {
    if (args.do_count == false) {
      if (args.format == csv)
        fprintf(args.output, "ID1,ID2,Distance\n");
      else
        fprintf(args.output, "Seq1,Seq2,Distance\n");
    }
  } else {
    distanceMatrix = new double[sequenceCount * sequenceCount];
    for (unsigned long i = 0; i < sequenceCount * sequenceCount; i++)
      distanceMatrix[i] = 100.;
    for (unsigned long i = 0; i < sequenceCount; i++)
      distanceMatrix[i * sequenceCount + i] = 0.;
  }

  long upperBound = (args.input2 && !do_fst) ? seqLengthInFile1 : sequenceCount,
       skipped_comparisons = 0;

  unsigned long global_hist[3][HISTOGRAM_BINS];
  for (unsigned long k = 0; k < 3; k++) {
    for (unsigned long r = 0; r < HISTOGRAM_BINS; r++) {
      global_hist[k][r] = 0;
    }
  }

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

  /*
   cout << recounts.length() << "/" << upperBound << endl;

  for (long k = 0; k < recounts.length(); k++) {
    cout << k << " = " << recounts.value(k) << endl;
  }
  */

  time_t before, after;
  time(&before);
  double weighted_counts[3] = {0., 0., 0.};
  /*
    if do_fst is true:
      index 0 -- file 1
      index 1 -- file 2
      index 2 -- file 1 vs file 2
  */

  bool cross_comparison_only = (args.input2 && !do_fst);

#pragma omp parallel shared(                                                   \
    skipped_comparisons, sequence_descriptors, resolutionOption, foundLinks,   \
    pairIndex, sequences, seqLengths, sequenceCount, firstSequenceLength,      \
    args, nameLengths, names, pairwise, percentDone, cerr, max, randFlag,      \
    distanceMatrix, upperBound, seqLengthInFile1, seqLengthInFile2, mean,      \
    randSeqs, weighted_counts, do_fst, randomized_fst, randomized_idx,         \
    recounts, cross_comparison_only, report_self) 

  {

    unsigned long histogram_counts[3][HISTOGRAM_BINS];
    for (unsigned long k = 0; k < 3; k++) {
      for (unsigned long r = 0; r < HISTOGRAM_BINS; r++) {
        histogram_counts[k][r] = 0;
      }
    }

#pragma omp for schedule(guided)
    for (long seq1 = 0; seq1 < upperBound; seq1++) {
      
      long mapped_id = randomized_fst ? randomized_idx.value(seq1)
                                      : (randSeqs ? randSeqs[seq1] : seq1);

      char *n1 = stringText(names, nameLengths, mapped_id),
           *s1 = stringText(sequences, seqLengths, mapped_id);

      long lowerBound = cross_comparison_only ? seqLengthInFile1 : seq1 + 1L,
           compsSkipped = 0, local_links_found = 0,
           instances1 =
               randomized_fst ? recounts.value(seq1) : counts.value(mapped_id);

      double local_max[3] = {0., 0., 0.}, local_sum[3] = {0., 0., 0.},
             local_weighted[3] = {0., 0., 0.};

      if (!cross_comparison_only) {
        if (instances1 > 1) {
          // add zero distances for self vs self
          unsigned long self_v_self = (instances1 * (instances1 - 1L)) >> 1;
          if (seq1 < seqLengthInFile1) {
            local_weighted[0] += self_v_self;
            histogram_counts[0][0] += self_v_self;
          } else {
            local_weighted[1] += self_v_self;
            histogram_counts[1][0] += self_v_self;
          }
        }
        if (report_self) {
          // if this is a single file run and -0 selected, then report distance
          // to self as 0
          if (args.format == csv) {
#pragma omp critical
            fprintf(args.output, "%s,%s,%g\n", n1, n1, 0.0);
          } else {
            if (args.format == csvn) {
#pragma omp critical
              fprintf(args.output, "%ld,%ld,%g\n", mapped_id, mapped_id, 0.0);
            }
          }
        }
      }

      for (unsigned long seq2 = lowerBound; seq2 < sequenceCount; seq2++) {
        long mapped_id2 = randomized_fst ? randomized_idx.value(seq2)
                                         : (randSeqs ? randSeqs[seq2] : seq2),
             which_bin = 0,
             weighted_count =
                 instances1 * (randomized_fst ? recounts.value(seq2)
                                              : counts.value(mapped_id2));

        if (do_fst) {
          if (seq1 < seqLengthInFile1) {
            if (seq2 >= seqLengthInFile1) {
              which_bin = 2; // file 1 vs file 2
            }
          } else {
            which_bin = 1;
          }
        }



        double thisD =
            sequence_descriptors
                ? computeTN93(s1, stringText(sequences, seqLengths, mapped_id2),
                              firstSequenceLength, resolutionOption, randFlag,
                              args.overlap, &(histogram_counts[which_bin][0]),
                              HISTOGRAM_SLICE, HISTOGRAM_BINS, weighted_count,
                              1L, &sequence_descriptors[mapped_id],
                              &sequence_descriptors[mapped_id2])
                : computeTN93(s1, stringText(sequences, seqLengths, mapped_id2),
                              firstSequenceLength, resolutionOption, randFlag,
                              args.overlap, &(histogram_counts[which_bin][0]),
                              HISTOGRAM_SLICE, HISTOGRAM_BINS, weighted_count);

       
        if (thisD >= -1.e-10 && thisD <= args.distance) {
          local_links_found += weighted_count;
          // char *s2 = stringText(sequences, seqLengths, seq1);
          if (!args.do_count) {
            if (args.format == csv) {
#pragma omp critical
              fprintf(args.output, "%s,%s,%g\n", n1,
                      stringText(names, nameLengths, mapped_id2), thisD);
            } else {
              if (args.format == csvn) {
#pragma omp critical
                fprintf(args.output, "%ld,%ld,%g\n", mapped_id, mapped_id2,
                        thisD);

              } else {
                distanceMatrix[mapped_id * sequenceCount + mapped_id2] = thisD;
                distanceMatrix[mapped_id2 * sequenceCount + mapped_id] = thisD;
              }
            }
          }
        }
        if (thisD <= -0.5) {
          compsSkipped += 1;
        } else {
          local_sum[which_bin] += thisD * (weighted_count);
          local_weighted[which_bin] += weighted_count;
          if (thisD > local_max[which_bin]) {
            local_max[which_bin] = thisD;
          }
        }
      }
#pragma omp critical
      {
        pairIndex += (args.input2 == NULL || do_fst)
                         ? (sequenceCount - seq1 - 1)
                         : seqLengthInFile2;
        foundLinks += local_links_found;
        skipped_comparisons += compsSkipped;
        for (int idx = 0; idx < 3; idx++) {
          if (local_max[idx] > max[idx]) {
            max[idx] = local_max[idx];
          }
          mean[idx] += local_sum[idx];
          weighted_counts[idx] += local_weighted[idx];
        }
      }

      if (!args.quiet && (pairIndex * 100. / pairwise - percentDone > 0.1 ||
                          seq1 == (long)sequenceCount - 1)) {
#pragma omp critical
        {
          time(&after);
          percentDone = pairIndex * 100. / pairwise;
          cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress"
                  ":"
               << setw(8) << percentDone << "% (" << setw(8) << foundLinks
               << " links found, " << setw(12) << std::setprecision(3)
               << pairIndex / difftime(after, before) << " evals/sec)";

          after = before;
        }
      }
    }
#pragma omp critical
    {
      for (unsigned long r = 0; r < 3; r++) {
        for (unsigned long k = 0; k < HISTOGRAM_BINS; k++) {
          global_hist[r][k] += histogram_counts[r][k];
        }
      }
    }
  } // end of the parallel loop

  if (args.do_count == false && args.format == hyphy) {
    fprintf(args.output, "{");
    for (unsigned long seq1 = 0; seq1 < sequenceCount; seq1++) {
      fprintf(args.output, "\n{%g", distanceMatrix[seq1 * sequenceCount]);
      for (unsigned long seq2 = 1; seq2 < sequenceCount; seq2++) {
        fprintf(args.output, ",%g",
                distanceMatrix[seq1 * sequenceCount + seq2]);
      }
      fprintf(args.output, "}");
    }
    fprintf(args.output, "\n}");

    delete[] distanceMatrix;
  }

  cerr << endl;

  ostream *outStream = &cout;

  if (args.output == stdout) {
    outStream = &cerr;
  }

  (*outStream) << "{" << endl
               << '\t' << "\"Actual comparisons performed\" :"
               << pairwise - skipped_comparisons << ',' << endl;
  (*outStream) << "\t\"Comparisons accounting for copy numbers \" :"
               << (weighted_counts[0] + weighted_counts[1] + weighted_counts[2])
               << ',' << endl;
  (*outStream) << "\t\"Total comparisons possible\" : " << pairwise << ','
               << endl;
  (*outStream) << "\t\"Links found\" : " << foundLinks << ',' << endl;
  (*outStream) << "\t\"Maximum distance\" : " << max[0] << ',' << endl;
  
  if (args.input2 == NULL) {
      (*outStream) << "\t\"Sequences\" : " << sequenceCount << ',' << endl;
  } else {
      (*outStream) << "\t\"Sequences in first file\" : " << seqLengthInFile1 << ',' << endl;
      (*outStream) << "\t\"Sequences in second file\" : " << seqLengthInFile2 << ',' << endl;
  }
  
  if (do_fst) {
    const char *keys[4] = {"File 1", "File 2", "Between", "Combined"};
    for (unsigned long k = 0; k < 3; k++) {
      (*outStream) << "\t\"Mean distance " << keys[k]
                   << "\" : " << mean[k] / weighted_counts[k] << ',' << endl;
    }
    double meta =
               (mean[0] + mean[1] + mean[2]) /
               (weighted_counts[0] + weighted_counts[1] + weighted_counts[2]),
           intra_pop =
               (mean[0] + mean[1]) / (weighted_counts[0] + weighted_counts[1]),
           pi_D = mean[2] / weighted_counts[2] - intra_pop;

    (*outStream) << "\t\"Mean distance [intra] "
                 << "\" : " << intra_pop << ',' << endl;
    (*outStream) << "\t\"Mean distance [combined]"
                 << "\" : " << meta << ',' << endl;
    (*outStream) << "\t\"F_ST\" : " << pi_D / (intra_pop + pi_D) << ',' << endl;

    for (unsigned long k = 0; k < 3; k++) {
      dump_histogram(outStream, keys[k], &(global_hist[k][0]));
      (*outStream) << ',' << endl;
    }
    for (unsigned long k = 0; k < HISTOGRAM_BINS; k++) {
      global_hist[0][k] += global_hist[1][k] + global_hist[2][k];
    }
    dump_histogram(outStream, keys[3], global_hist[0]);
  } else {
    (*outStream) << "\t\"Mean distance\" : " << mean[0] / weighted_counts[0]
                 << ',' << endl;
    dump_histogram(outStream, NULL, global_hist[0]);
  }
  (*outStream) << '}' << endl;

  if (args.do_count) {
    fprintf(args.output, "Found %ld links among %ld pairwise comparisons\n",
            foundLinks, pairwise - skipped_comparisons);
  }

  if (randFlag)
    delete[] randFlag;

  if (randSeqs)
    delete[] randSeqs;
  else
    delete[] sequence_descriptors;

  return 0;
}
