
#include "tn93_shared.h"
#include "argparse_cf.hpp"

using namespace std;
using namespace argparse;

//---------------------------------------------------------------

void handle_a_sequence (StringBuffer& current_sequence, double * resolutions, const long firstSequenceLength, long instance_count, bool ignore_ambigs) {
    long resolution_index = 0;
    const unsigned char *s = (const unsigned char *)current_sequence.getString();
    double local_instance_count = instance_count;
    for (long c = 0; c < firstSequenceLength; c++, resolution_index+=4) {
        const long * res = resolve_char(s[c], false, ignore_ambigs);
        if (ignore_ambigs == false) {
          local_instance_count = instance_count * resolution_count(s[c], false);
        }
        resolutions [resolution_index]   += local_instance_count*res[0];
        resolutions [resolution_index+1] += local_instance_count*res[1];
        resolutions [resolution_index+2] += local_instance_count*res[2];
        resolutions [resolution_index+3] += local_instance_count*res[3];
    }
}   

//---------------------------------------------------------------

void handle_a_sequence_aa (StringBuffer& current_sequence, double * resolutions, const long firstSequenceLength, long instance_count, bool ignore_ambigs) {
    long resolution_index = 0;
    const unsigned char *s = (const unsigned char *)current_sequence.getString();
    double local_instance_count = instance_count;
    for (long c = 0; c < firstSequenceLength; c++, resolution_index+=20) {
        const long * res = resolve_char(s[c], true, ignore_ambigs);
        if (ignore_ambigs == false) {
          local_instance_count = instance_count * resolution_count(s[c], true);
        }
        for (long k = 0; k < 20; k++) {
          resolutions [resolution_index + k]   += local_instance_count*res[k];
        }
    }
}   


//---------------------------------------------------------------

void exportJSON (FILE* out, long rows, long columns, const double* counts, bool do_prot) {
  fprintf (out, "{");
  
  if (counts) {
    bool do_comma = false;
    for (long r = 0; r < rows; r++) {
      bool first = true;
      for (long c = 0; c < columns; c++) {
        if (counts [r*columns + c] > 0.0) {
          if (first) {
            fprintf (out, "%c\n\t\"%ld\" : {\"%c\" : %g", do_comma ? ',' : ' ', r + 1, unmap_char(c,do_prot), counts[r*columns+c]);
            do_comma = true;
            first = false;
          } else {
            fprintf (out, ", \"%c\" : %g", unmap_char(c,do_prot) , counts[r*columns+c]);
          }
        }
      }
      if (!first) {
        fprintf (out, "}");
      }
    }
  }
  
  fprintf (out, "\n}\n");
}

//---------------------------------------------------------------

int main (int argc, const char * argv[]) {

    args_t args = args_t (argc, argv);
    initAlphabets(args.data == protein);

    char automatonState = 0,
         fasta_result = 2;
    // 0 - between sequences
    // 1 - reading sequence name
    // 2 - reading sequence
    
    double * freq_counts = NULL;
    Vector copy_count;
    copy_count.appendValue(0);

 
    long sequences_read = 0,
         firstSequenceLength = 0,
         dim = (args.data == dna) ? 4: 20;
          
    while (fasta_result == 2) {
        
        fasta_result = readFASTA (args.input, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength, true, &copy_count, args.counts_in_name);
        if (fasta_result == 1) {
            return 1;
        }
        if (sequences_read == 0) {
            freq_counts = new double [dim*firstSequenceLength];
            for (long k = 0; k < dim*firstSequenceLength; k++) freq_counts [k] = 0.;
        }
        
        if (args.data == protein) 
          handle_a_sequence_aa (sequences, freq_counts, firstSequenceLength, copy_count.value(0), args.ambig == ignore);        
        else
          handle_a_sequence (sequences, freq_counts, firstSequenceLength, copy_count.value(0), args.ambig == ignore);
        sequences_read ++;
        if (args.quiet == false && sequences_read % 1000 == 0) {
            cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << sequences_read << " reads";
        }

    }
    
    if (args.quiet == false) {
      cerr << endl;
    }
    
    exportJSON (args.output, firstSequenceLength, dim, freq_counts, args.data == protein);
    if (freq_counts) delete [] freq_counts;
    
    return 0;

}

