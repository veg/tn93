
#include "tn93_shared.h"
#include "argparse_fasta_diff.hpp"
#include <map>
#include <list>
#include <string>

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

int main(int argc, const char *argv[]) {
    args_t args = args_t(argc, argv);
    
    
    long firstSequenceLength = 1L;
    
    try {
        char automatonState = 0;
        // 0 - between sequences
        // 1 - reading sequence name
        // 2 - reading sequence
        
        std::map<std::string, std::string> sequences_to_add;
        
        /*auto compare_records = [&] (const std::string & id1, const char* seq_value) -> int {
            std::map<std::string, std::string>::iterator it = sequences_to_add.find (id1);
            if (it !=  sequences_to_add.end()) {
                if (args.checks == id) {
                    return 1;
                }
                return it->second == std::string (seq_value) ? 2 : 3;
            }
            return 0;
        };*/
        


        initAlphabets(false, NULL, true);
        long non_unique_ids = 0;
        bool sequence_deduplication = false;
      
      if (args.checks == argparse::sequence) {
        if (args.op == argparse::remove) {
          sequence_deduplication = true;
        } else {
          throw (std::string ("sequence matching mode MUST also use the remove operation"));
        }
      }
      

        while (long state = readFASTA(args.input_add, automatonState, names, sequences, nameLengths,
                         seqLengths, firstSequenceLength, true)) { // read sequences one by one
            
            if (state == 2 || state == 3) {
                if (sequence_deduplication) {
                  // for this run mode, the map is inverted : sequence  -->  ID
                  sequences_to_add[sequences.getString()] = names.getString();
                } else {
                  long current_size = sequences_to_add.size();
                  std::string check_existing = sequences_to_add[names.getString()];
                  sequences_to_add[names.getString()] = sequences.getString();
                  if (current_size == sequences_to_add.size()) {
                      if (check_existing == sequences_to_add[names.getString()]) {
                          //cerr << std::string (names.getString()) + " is not a unique sequence ID in the master file" << endl;
                          non_unique_ids++;
                      } else {
                          throw (std::string (names.getString()) + " is not a unique sequence ID in the master file and different instances of this name have different sequences");
                      }
                  }
                }
                
                if (state == 3) {
                    break;
                }
            } else {
                throw (std::string ("Error reading FASTA file"));
            }
            //dump_fasta (sequences.getString(), firstSequenceLength, stdout);
        }
        
        long sequences_to_add_count = sequences_to_add.size(),
             duplicates = 0L,
             master_sequences = 0L,
             output = 0L;

        auto echo_fasta_sequence = [&] (const char* id, const char* data, FILE * where) -> void {
            output ++;
            fprintf (where, ">%s\n%s\n", id, data);
        };


        /*
         scan master records one by one;
         output those that are not present in sequences_to_add
         those that are present in sequences_to_add will be
        */
        
        std::list <std::string> updated_sequences;
        
        while (long state = readFASTA(args.input_master, automatonState, names, sequences, nameLengths, seqLengths, firstSequenceLength, true)) {
            if (state == 2 || state == 3) {
                std::string master_id (names.getString());
                //cerr << master_id << endl;
                master_sequences++;
              
                bool records_different,
                     ids_different;
                if (!sequence_deduplication) {
                  auto it = sequences_to_add.find (master_id);
                  ids_different = it == sequences_to_add.end();
                  
                  if (ids_different) {
                      records_different = true;
                  } else {
                      if (args.checks == argparse::id_and_sequence) {
                          records_different = it->second != std::string (sequences.getString());
                      } else {
                          records_different = false;
                      }
                  }
                } else {
                  ids_different = true;
                  auto it = sequences_to_add.find (sequences.getString());
                  records_different = it == sequences_to_add.end();
                }
                
                //cout << records_different << " " << ids_different << endl;
                
                if (!records_different) { // not unique
                    duplicates++;
                    if (args.op == argparse::replace) {
                        if (state == 3) {
                            break;
                        }
                        continue;
                    }
                    if (args.op == argparse::add) {
                        echo_fasta_sequence (master_id.c_str(), sequences.getString(), args.output);
                    }
                    sequences_to_add.erase (master_id);
                } else {
                    if (ids_different || args.op != argparse::replace) {
                        if (!ids_different) {
                            updated_sequences.push_back (master_id);
                            if (args.op == argparse::add) {
                                sequences_to_add.erase (master_id);
                            }
                        }
                        echo_fasta_sequence (master_id.c_str(), sequences.getString(), args.output);
                    }
                }
                
                if (state == 3) {
                    break;
                }
            } else {
                throw (std::string ("Error reading FASTA file"));
            }
        }
                
                // now output remaining sequences from the FASTA file
        
        if (args.op != argparse::remove) {
            for (std::pair<std::string, std::string> it : sequences_to_add) {
                echo_fasta_sequence (it.first.c_str(), it.second.c_str(), args.output);
            }
        }
        funlockfile (args.output);

        if (!args.quiet) {
            cerr << endl << "{" << endl << "'master' :" << master_sequences << ',' << endl <<
                      "'fasta'  : " << sequences_to_add_count << ',' << endl <<
                      "'output' : " << output << ',' << endl <<
                      "'non-unqiue' : " << non_unique_ids << ',' << endl <<
                      "'duplicate' : " << duplicates << ',' << endl <<
                      "'updated_sequences' : [";
            
            long k = 0;
            for (std::string s : updated_sequences) {
                if (k) {
                    cerr << ", ";
                }
                cerr << "'" << s << "'";
                k++;
            }
            cerr << "]" << endl << "}" << endl;
        }
        
        /*
        if (!args.quiet) {
            cerr << "Read " << sequences_to_add.size() << " sequences from the FASTA file" << endl;
        }*/

    } catch (const std::string err) {
        cerr << "ERROR: " << err << endl;
        return 1;
    }
    
    return 0;
}
