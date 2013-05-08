
#include "tn93_shared.h"
#include <float.h>
#include <math.h>

using namespace std;

static char Usage[] = "nucfreqsfasta"
                      "\n\t<FASTA file OR - for stdin >"
                      "\n\t<output file OR - for stdout>"
                      "\n\tNORMAL|BINMIX -- output format choices"
                      "\n\t[background rate, [p-value]]"
                        ;

;

double computeBinomialProb (double p, long coverage, long mutation) {
    double sum = 0.0,
           term_update = p/(1.-p),
           cov = coverage,
           current_term = pow ((1.-p),cov),
           idx = 0.;
    
    long   index = 0;
    
    while  (current_term > 0. && index < mutation) {
        sum +=current_term;
        current_term *= term_update * (cov-idx) / (idx+1);
        index ++;
        idx += 1.;
    }

    return 1.-sum;
}

//---------------------------------------------------------------

void handle_a_sequence (StringBuffer& current_sequence, long * resolutions, const long firstSequenceLength) {
    long resolution_index = 0;
    const unsigned char *s = (const unsigned char *)current_sequence.getString();
    for (long c = 0; c < firstSequenceLength; c++, resolution_index+=4) {
        const long * res = resolve_char(s[c]);
        resolutions [resolution_index] += res[0];
        resolutions [resolution_index+1] += res[1];
        resolutions [resolution_index+2] += res[2];
        resolutions [resolution_index+3] += res[3];
    }
}   



//---------------------------------------------------------------

int main (int argc, const char * argv[])
{
    if (argc != 4 && argc != 5 && argc != 6)
    {
        cerr << "Usage is `" << Usage << "'." << endl;
        exit(1);
    }

    const    char *S = argv[1];
    bool     bin_mix = strcmp (argv[3], "BINMIX") == 0;
    
    double   filter_rate = -1.,
             p_value = 0.01;

    long     firstSequenceLength  = 0,
            * freq_counts = NULL;
    
    if (argc >=5 ) {
        filter_rate = atof (argv[4]);
        if (filter_rate <= 0. || filter_rate >= 1.) {
            cerr << "Invalid background rate supplied [must be in (0,1)]: `" << argv[4] << "'." << endl;
            return 1;
           
        }
        if (argc == 6) {
            p_value = atof (argv[5]);
            if (filter_rate <= 0. || filter_rate >= 1.) {
                cerr << "Invalid p-value supplied [must be in (0,1)]: `" << argv[5] << "'." << endl;
                return 1;
                
            }
        }
    }
    
    FILE *F  = strcmp (S,"-") == 0 ? stdin: fopen(S, "r"),
          *FO = strcmp (argv[2],"-") == 0? stdout: fopen (argv[2], "w");

    if (F == NULL) {
        cerr << "Cannot open file `" << S << "'." << endl;
        return 1;
    }

    if (FO == NULL) {
        cerr << "Cannot open file `" << FO << "'." << endl;
        return 1;
    }
    
    //cerr << computeBinomialProb(0.0025, 100, 3) << endl;

 
    char automatonState = 0,
         fasta_result = 2;
    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence

    initAlphabets();
 
    long sequences_read = 0;
          
    while (fasta_result == 2) {
        fasta_result = readFASTA (F, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength, true);
        if (fasta_result == 1) {
            return 1;
        }
        if (sequences_read == 0) {
            freq_counts = new long [4*firstSequenceLength];
            for (long k = 0; k < 4*firstSequenceLength; k++) freq_counts [k] = 0;
        }
        
        handle_a_sequence (sequences, freq_counts, firstSequenceLength);
        sequences_read ++;
        if (sequences_read % 1000 == 0) {
            cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << sequences_read << " reads";
        }

    }
    handle_a_sequence (sequences, freq_counts, firstSequenceLength);
    cerr << endl;
    fclose(F);
    
    if (freq_counts) {
        delete [] freq_counts;
    }
    
    if (bin_mix) {
         for (long s = 0; s < firstSequenceLength; s++){
            long coverage = freq_counts[s*4] + freq_counts[s*4+1] + freq_counts[s*4+2] + freq_counts[s*4+3],
            conensus     = 0,
            consensus_id = 0;
            
            for (long i = 0; i < 4; i++) {
                if (freq_counts[s*4+i] > conensus) {
                    conensus = freq_counts[s*4+i];
                    consensus_id = i;
                }
            }

            for (long i = 0; i < 4; i++) {
                if (i!=consensus_id && freq_counts[s*4+i]) {
                    fprintf (FO, "%ld\t%ld\n", coverage, coverage-freq_counts[s*4+i]);
                }
            }
        }
    } else {
        if (filter_rate < 0.0){
            for (long s = 0; s < firstSequenceLength; s++){
                long coverage = freq_counts[s*4] + freq_counts[s*4+1] + freq_counts[s*4+2] + freq_counts[s*4+3];
                fprintf (FO, "%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",s+1,coverage,freq_counts[s*4],freq_counts[s*4+1],freq_counts[s*4+2] ,freq_counts[s*4+3]);
            }
           
        } else {
            char nucs [5] = "ACGT";
            for (long s = 0; s < firstSequenceLength; s++){
                long coverage = freq_counts[s*4] + freq_counts[s*4+1] + freq_counts[s*4+2] + freq_counts[s*4+3];
                
                for (long k = 0; k < 4; k++) {
                    double supp = computeBinomialProb (filter_rate, coverage, freq_counts[s*4+k]);
                    if (supp < p_value) {
                        fprintf (FO, "%ld\t%ld\t%ld\t%c:%g\n",s+1,coverage,freq_counts[s*4+k], nucs[k],supp);
                    }
                }
             }
        }
    }

    if (FO != stdout)
        fclose (FO);
    return 0;

}

