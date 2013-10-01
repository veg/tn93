
#include "tn93_shared.h"

#ifdef  _OPENMP
  #include "omp.h"
#endif

#include "argparse.hpp"


using namespace std;

#define HISTOGRAM_BINS 200
#define HISTOGRAM_SLICE ((double)HISTOGRAM_BINS)

static char Usage[] = "TN93dist"
                      "\n\t<FASTA file OR - for stdin >"
                      "\n\t<[output file OR - for stdout> OR COUNT>"
                      "\n\t<distance threshold>"
                      "\n\t<how to handle ambiguities; one of RESOLVE, AVERAGE, SKIP, GAPMM>"
                      "\n\t<output format; one of CSV, CSVN (numeric IDs instead of sequence names) HYPHY>"
                      "\n\t<minimum overlap between sequences: integer >= 1>"
                      "\n\t[BOOTSTRAP 0 or 1]"
                      "\n\t[SECOND FILE]";



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


    FILE *F  = strcmp (S,"-") == 0 ? stdin: fopen(S, "r"),
          *FO = count_only ? NULL : (strcmp (argv[2],"-") == 0? stdout: fopen (argv[2], "w")),
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

    if (strcmp (argv[4], "GAPMM") == 0) {
        resolutionOption = GAPMM;
    } else if (strcmp (argv[4], "AVERAGE") == 0) {
        resolutionOption = AVERAGE;
    } else if (strcmp (argv[4], "SKIP") == 0) {
        resolutionOption = SKIP;
    }

    StringBuffer names,
                 sequences;

    Vector       nameLengths,
                 seqLengths,
                 counts;



    char automatonState = 0;
    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence

    nameLengths.appendValue (0);
    seqLengths.appendValue (0);
    initAlphabets();
    if (readFASTA (F, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength) == 1)
        return 1;

    fclose(F);

    unsigned long seqLengthInFile1 = seqLengths.length()-1,
                  seqLengthInFile2 = 0;

    if (F2) {
        if (readFASTA (F2, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength) == 1)
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
           max         = 0.0,
           mean        = 0.0;


    unsigned char    doCSV = 1;

    if (strcmp (argv[5], "HYPHY") == 0)
        doCSV = 0;
    else if (strcmp (argv[5], "CSVN") == 0)
        doCSV = 2;


    long * randFlag = NULL;
    long * randSeqs = NULL;

    if (argc >= 8 && atoi (argv[7]) > 0) {
        if (argc == 8) {
            cerr << endl << "Randomizing site order..." << endl;
            randFlag = new long [firstSequenceLength];
            bool   *included = new bool [firstSequenceLength];
            for (long k = 0; k < firstSequenceLength; k++) {
                included[k] = false;
            }
            init_genrand (time(NULL) + getpid ());

            unsigned long normalizer = RAND_RANGE / firstSequenceLength,
                          max_site   = 0;


            for (unsigned long i = 0; i < firstSequenceLength; i++) {
                randFlag[i] = genrand_int32 () / normalizer;
                cerr << randFlag[i] << " ";
                included[randFlag[i]] = true;
                if (randFlag[i] > max_site) max_site = randFlag[i];
            }

            long total_resampled = 0;
            for (long k = 0; k < firstSequenceLength; k++) {
                total_resampled += included[k];
            }

            delete [] included;
            cerr << endl << "Unique sites included in the resampled order " << total_resampled << endl;
        } else {
            if (argc == 9) {
                randSeqs = new long [sequenceCount];

                for (unsigned long k = 0; k < sequenceCount; k++) {
                    randSeqs [k] = k;
                }

                init_genrand (time(NULL) + getpid ());


                for (unsigned long i = 0; i < sequenceCount; i++) {
                    long id = genrand_int32 () / (RAND_RANGE / (sequenceCount-i)),
                         t = randSeqs[i+id];
                    randSeqs[i+id]= randSeqs[i];
                    randSeqs[i] = t;
                }



                cerr << endl << "Randomized sequence assignment to datasets " << endl;
            }
        }
    }



    cerr << endl << "Progress: " << setw(8) << 0.0 << "% (" << setw(8) << 0 << " links found)";

    long pairIndex  = 0,
         foundLinks = 0;

    double *distanceMatrix = NULL;

    if (doCSV) {
        if (!count_only) {
            if (doCSV > 1)
                fprintf (FO, "ID1,ID2,Distance\n");
            else
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
         
    
    unsigned long global_hist [HISTOGRAM_BINS];
    for (unsigned long k = 0; k < HISTOGRAM_BINS; k++) {
      global_hist[k] = 0;
    }
    

    #pragma omp parallel shared(count_only, skipped_comparisons, resolutionOption, foundLinks,pairIndex,sequences,seqLengths,sequenceCount,firstSequenceLength,distanceThreshold, nameLengths, names, pairwise, percentDone,FO,cerr,max,randFlag,doCSV,distanceMatrix, upperBound, argc, seqLengthInFile1, seqLengthInFile2, min_overlap, mean, randSeqs)
    
    {
    
      unsigned long histogram_counts [HISTOGRAM_BINS];
      for (unsigned long k = 0; k < HISTOGRAM_BINS; k++) {
        histogram_counts[k] = 0;
      }
    
      #pragma omp for schedule (dynamic)
      for (long seq1 = 0; seq1 < upperBound; seq1 ++)
      {
          long mapped_id = randSeqs?randSeqs[seq1]:seq1;

          char *n1 = stringText (names, nameLengths, mapped_id),
                *s1 = stringText(sequences, seqLengths, mapped_id);

          long lowerBound = argc == 9 ? seqLengthInFile1 : seq1 +1,
               compsSkipped      = 0,
               local_links_found = 0;

          double local_max = 0.,
                 local_sum = 0.;

          for (unsigned long seq2 = lowerBound; seq2 < sequenceCount; seq2 ++)
          {
              long mapped_id2 = randSeqs?randSeqs[seq2]:seq2;
              double thisD = computeTN93(s1, stringText(sequences, seqLengths, mapped_id2), firstSequenceLength, resolutionOption, randFlag, min_overlap, histogram_counts, HISTOGRAM_SLICE, HISTOGRAM_BINS);

              if (thisD >= -1e-10 && thisD <= distanceThreshold)
              {
                  local_links_found ++;
                  //char *s2 = stringText(sequences, seqLengths, seq1);
                  if (!count_only) {
                      if (doCSV == 1) {
                          #pragma omp critical
                          fprintf (FO,"%s,%s,%g\n", n1, stringText (names, nameLengths, mapped_id2), thisD);
                      } else {
                          if (doCSV == 2) {
                              #pragma omp critical
                              fprintf (FO,"%ld,%ld,%g\n", mapped_id, mapped_id2, thisD);

                          } else {
                              distanceMatrix[mapped_id*sequenceCount+mapped_id2] = thisD;
                              distanceMatrix[mapped_id2*sequenceCount+mapped_id] = thisD;
                          }
                      }
                  }
              }
              if (thisD <= -0.5) {
                  compsSkipped ++;
              }
              else {
                  local_sum += thisD;
                  if (thisD > local_max) {
                      local_max = thisD;
                  }
              }


          }
          #pragma omp critical
          {
              foundLinks += local_links_found;
              skipped_comparisons += compsSkipped;
              if (local_max > max) {
                  max = local_max;
              }
              pairIndex += (argc < 9) ? (sequenceCount - seq1 - 1) : seqLengthInFile2;
              mean += local_sum;
          }

          if (pairIndex * 100. / pairwise - percentDone > 0.1 || seq1 == (long)sequenceCount - 1)
          {
              #pragma omp critical
              percentDone = pairIndex * 100. / pairwise;
              #pragma omp critical
              cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << percentDone << "% (" << setw(8) << foundLinks << " links found)";
          }
      }
      #pragma omp critical
      {
        for (unsigned long k = 0; k < HISTOGRAM_BINS; k++) {
          global_hist[k] += histogram_counts[k];
        }
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
    
    ostream * outStream = &cout;
    
    if (FO == stdout) {
      outStream = &cerr;
    }
    
    (*outStream) << "{" << endl << '\t' << "\"Actual comparisons performed\" :" << pairwise-skipped_comparisons << ',' << endl;
    (*outStream) << "\t\"Total comparisons possible\" : " << pairwise << ',' << endl;
    (*outStream) << "\t\"Links found\" : " << foundLinks << ',' << endl;
    (*outStream) << "\t\"Maximum distance\" : " << max << ',' << endl;
    (*outStream) << "\t\"Mean distance\" : " << mean/(pairwise-skipped_comparisons) << ',' << endl;
    (*outStream) << "\t\"Histogram\" : [";
    for (unsigned long k = 0; k < HISTOGRAM_BINS; k++) {
      if (k) {
        (*outStream) << ',';
      }
      (*outStream) << '[' << (k+1.) / HISTOGRAM_BINS << ',' << global_hist[k] << ']';
    }
    (*outStream) << "]" << endl;
    
    (*outStream) << '}' << endl;
    
    if (count_only) {
        cout << "Found " << foundLinks << " links among " << pairwise-skipped_comparisons << " pairwise comparisons" << endl;
    }

    if (randFlag)
        delete [] randFlag;

    if (randSeqs)
        delete [] randSeqs;

    if (FO)
        fclose (FO);
    return 0;

}

