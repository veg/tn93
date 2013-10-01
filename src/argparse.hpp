
#ifndef ARGPARSE_H
#define ARGPARSE_H

// argument defaults

#define PROGNAME "tn93"
#define DEFAULT_AMBIG             RESOLVE
#define DEFAULT_FORMAT            CSV
#define DEFAULT_DISTANCE          0.015
#define DEFAULT_COUNTS_IN_NAME    ':'
#define DEFAULT_OVERLAP           100UL

namespace argparse
{
    enum ambig_t {
        RESOLVE,
        AVERAGE,
        SKIP,
        GAPMM
    };
    
    enum format_t {
      CSV,
      CSV_N,
      HYPHY
    };

    class args_t
    {
    public:
 
        FILE            * output,
                        * input1,
                        * input2;
             
        double          distance;
        ambig_t         ambig;
        format_t        format;
        unsigned long   overlap;
        bool            do_bootstrap;
        bool            do_count;
        char            counts_in_name;
        
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input    ( const char * );
        void parse_second_in( const char * );
        void parse_output   ( const char * );
        void parse_distance ( const char * );
        void parse_overlap  ( const char * );
        void parse_format   ( const char * );
        void parse_counts_in_name   ( const char * );
        void parse_bootstrap( void );
        
    };
}

#endif // ARGPARSE_H
