
#ifndef ARGPARSE_H
#define ARGPARSE_H
#include <cstdarg>

// argument defaults


#define PROGNAME                 "fasta_diff"

#define        DEFAULT_OPERATION         add
#define        DEFAULT_MATCH             id


namespace argparse {
    void ERROR( const char * msg, ... );
    
    enum file_operation {
        add,
        replace,
        remove
    };
    
    enum match_mode {
        id,
        id_and_sequence,
        sequence
    };
    
    
    class args_t {
    public:
        
        FILE * input_master,
        * input_add,
        * output;
        
        file_operation      op;
        match_mode          checks;
        
        bool                quiet;
        
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input_master       ( const char * );
        void parse_input_add          ( const char * );
        void parse_output             ( const char * );
        void parse_quiet              ( void );
        void parse_file_operation     ( const char * );
        void parse_match_mode         ( const char * );
    };
}

#endif // ARGPARSE_H
