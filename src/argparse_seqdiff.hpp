
#ifndef ARGPARSE_H
#define ARGPARSE_H
#include <cstdarg>

// argument defaults


#define PROGNAME                 "seq_diff"

#define        DEFAULT_AMBIG             informative
#define        DEFULT_FORMAT             json


namespace argparse {
    void ERROR( const char * msg, ... );
    
    enum output_format {
        json,
        csv
    };
    
    enum ambig_mode {
        all,
        informative,
        resolve
    };
    
    
    class args_t {
    public:
        
        FILE * input,
        * reference,
        * output;
        
        output_format       out_format;
        ambig_mode          ambig;
        
        bool                quiet;
        bool                compress;
        
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input              ( const char * );
        void parse_output             ( const char * );
        void parse_quiet              ( void );
        void parse_compress           ( void );
        void parse_reference          ( const char * );
        void parse_format             ( const char * );
        void parse_ambig              ( const char * );
    };
}

#endif // ARGPARSE_H
