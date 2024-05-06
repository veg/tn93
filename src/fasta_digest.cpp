
#include "tn93_shared.h"
#include "argparse_fasta_digest.hpp"
#include <map>
#include <list>
#include <string>

using namespace std;
using namespace argparse;

Vector* patterns[13] = {NULL};

// bit masks for
const unsigned long masks [16] = {
    0b00000000000000000000000000000011,
    0b00000000000000000000000000001100,
    0b00000000000000000000000000110000,
    0b00000000000000000000000011000000,
    0b00000000000000000000001100000000,
    0b00000000000000000000110000000000,
    0b00000000000000000011000000000000,
    0b00000000000000001100000000000000,
    0b00000000000000110000000000000000,
    0b00000000000011000000000000000000,
    0b00000000001100000000000000000000,
    0b00000000110000000000000000000000,
    0b00000011000000000000000000000000,
    0b00001100000000000000000000000000,
    0b00110000000000000000000000000000,
    0b11000000000000000000000000000000
};

int  allowed_chars [255];

unsigned patternToCode (const StringBuffer& pattern) {
    
    unsigned code = 0,
             L = pattern.length();
    
    for (int i = 0; i < L; i++) {
        code += allowed_chars[pattern.getChar (i)] << (2*i);
    }
    
    cout << pattern.getString() << "=" << code << endl;
    return code;
}

//---------------------------------------------------------------

int main(int argc, const char *argv[]) {
    args_t args = args_t(argc, argv);
    
    initAlphabets(false, NULL, true);
    
    
    for (int i = 0; i < 256; i++) {
        allowed_chars[i] = -1;
    }
    
    allowed_chars['A'] = 0;
    allowed_chars['C'] = 1;
    allowed_chars['G'] = 2;
    allowed_chars['T'] = 3;
    
    char complement_chars[255] = {0};
    complement_chars['A'] = 'T';
    complement_chars['T'] = 'A';
    complement_chars['C'] = 'G';
    complement_chars['G'] = 'C';


    auto process_motif = [&] (StringBuffer& pattern) -> void {
        int L = pattern.length();
        if (L >= 4 && L <= 16) {
            cout << pattern.getString () << endl;
            if (!patterns[L]) {
                patterns[L] = new Vector;
            }
            patterns[L]->appendValue(patternToCode (pattern));
            
            if (args.rc == complement) {
                StringBuffer rc;
                for (int i = L-1; i>=0; i--) {
                    rc.appendChar (complement_chars[pattern.getChar (i)]);
                }
                patterns[L]->appendValue(patternToCode (rc));
            }
            
            pattern.resetString ();
        } else {
            throw (std::string ("Motif '") + std::string (pattern.getString()) + "' is not 4-16 nucleotides in length");
        }
    };
    
    try {
        // parse the patterns
        StringBuffer pattern;
        char automatonState = 0;
        long firstSequenceLength = strlen (args.motif_list);
 
        for (unsigned i = 0; i < firstSequenceLength; i++) {
            char current_char = toupper (args.motif_list[i]);
            if (allowed_chars [current_char] >= 0) {
                automatonState = 1;
                pattern.appendChar (current_char);
            } else {
                if (automatonState == 1) {
                    process_motif (pattern);
                }
                automatonState = 0;
            }
        }
        
        if (automatonState == 1) {
            process_motif (pattern);
        }

        return 1;
        automatonState = 0;
        while (long state = readFASTA(args.input, automatonState, names, sequences, nameLengths,
                                      seqLengths, firstSequenceLength, true)) { // read sequences one by one
            
            firstSequenceLength = 0L;
        }

    } catch (const std::string err) {
        cerr << "ERROR: " << err << endl;
        return 1;
    }
    
    return 0;
}
