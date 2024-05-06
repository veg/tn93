
#ifndef ARGPARSE_H
#define ARGPARSE_H
#include <cstdarg>


// argument defaults


#define PROGNAME                         "fasta_digest"

#define        DEFAULT_RC                complement
#define        DEFAULT_MATCH             exact


namespace argparse {
    void ERROR( const char * msg, ... );
    
    enum rc_option {
        complement,
        nocomplement
    };
    
    enum match_mode {
        exact,
        oneoff,
        transition
    };
    
    
    class args_t {
    public:
        
        FILE * input,
             * output;
        
        char * motif_list;
        
        rc_option           rc;
        match_mode          match;
        
        
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input              ( const char * );
        void parse_enzyme_lis         ( const char * );
        void parse_output             ( const char * );
        void parse_rc_option          ( const char * );
        void parse_match_mode         ( const char * );
        void parse_motif              ( const char * );
   };
}

#endif // ARGPARSE_H
