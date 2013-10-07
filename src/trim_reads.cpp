
#include "tn93_shared.h"
#include "argparse_trim.hpp"

using namespace std;
using namespace argparse;

//---------------------------------------------------------------

void handle_a_sequence (FILE* output, StringBuffer& current_name, StringBuffer& current_sequence, unsigned long from, unsigned long to, double coverage, ambig_mode ambig_mode, bool is_protein) {
    long covered = 0;
    const unsigned char *s = (const unsigned char *)current_sequence.getString();
    for (unsigned long c = from; c <= to; c++) {
        double rc = resolution_count (s[c], false);
        if (rc == 1.) {
          covered ++;
        } else if (rc > 0.) {
          switch (ambig_mode) {
            case nfold: 
              if (rc > (is_protein? 0.05 : 0.25)) {
                covered ++;
              }
              break;
            case threefold:
              if (rc > 1./3.) {
                covered ++;
              }
              break;
            case any:
            case gaponly:
              break;
          }
        }
    }
    if (covered/ (to-from+1.0) >= coverage) {
      dump_sequence_fasta (0, output, 0, NULL,  is_protein, from, to);
    }
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
    
    Vector copy_count;
    copy_count.appendValue(0);

 
    long sequences_read = 0,
         firstSequenceLength = 0;
          
    while (fasta_result == 2) {
        
        fasta_result = readFASTA (args.input, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength, true);
        if (fasta_result == 1) {
            return 1;
        }
        if (sequences_read == 0) {
              if (args.start_coord >= firstSequenceLength) {
                ERROR ("start_coord must be less than the sequence length (%ld), had (%ld)", firstSequenceLength, args.start_coord);
              }
              if (args.end_coord >= firstSequenceLength) {
                args.end_coord = firstSequenceLength-1;
              }
              if (args.end_coord <= args.start_coord) {
                ERROR( "start of the filtering frame must be less than end of the frame, had: %ld - %ld", args.start_coord, args.end_coord );
              }
        }
        
  
        handle_a_sequence (args.output, names, sequences, args.start_coord, args.end_coord, args.coverage, args.ambig, args.data == protein);
        
        sequences_read ++;
        if (args.quiet == false && sequences_read % 1000 == 0) {
            cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(8) << sequences_read << " reads";
        }

    }
    
    if (args.quiet == false) {
      cerr << endl;
    }
    
    
    return 0;

}

