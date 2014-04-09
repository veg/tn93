
#ifndef ARGPARSE_H
#define ARGPARSE_H

// argument defaults


#define PROGNAME                 "tn93"
#define DEFAULT_AMBIG             resolve
#define DEFAULT_FORMAT            csv
#define DEFAULT_DISTANCE          0.015
#define DEFAULT_COUNTS_IN_NAME    ':'
#define DEFAULT_OVERLAP           100
#define DEFAULT_INCLUDE_PROB      1.0


namespace argparse
{
  
    enum format_t {
      csv,
      csvn,
      hyphy
    };

    enum ambig_t {
      resolve,
      average,
      skip,
      gapmm,
      subset
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
        bool            quiet;
        bool            do_fst;
        char            counts_in_name;
        double          include_prob;
        char            *ambigs_to_resolve;
      
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input    ( const char * );
        void parse_second_in( const char * );
        void parse_output   ( const char * );
        void parse_distance ( const char * );
        void parse_overlap  ( const char * );
        void parse_format   ( const char * );
        void parse_ambig    ( const char * );
        void parse_include_prob   ( const char * );
        void parse_counts_in_name
                            ( const char * );
        void parse_bootstrap( void );
        void parse_count    ( void );
        void parse_quiet    ( void );
        void parse_fst      ( void );
      
    };
}

#endif // ARGPARSE_H
