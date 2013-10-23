
#ifndef ARGPARSE_H
#define ARGPARSE_H

// argument defaults


#define PROGNAME                 "seqcoverage"
#define DEFAULT_AMBIG             ignore
#define DEFAULT_DATA              dna
#define DEFAULT_COUNTS_IN_NAME    ':'

namespace argparse
{
  
    enum ambig_t {
      ignore,
      average,
    };
    
    enum data_t {
      dna,
      protein
    };

  class args_t
    {
    public:
 
        FILE            * output,
                        * input;
             
        ambig_t         ambig;
        data_t          data;
        char            counts_in_name;
        bool            quiet;
        
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input    ( const char * );
        void parse_output   ( const char * );
        void parse_ambig    ( const char * );
        void parse_data     ( const char * );      
        void parse_counts_in_name ( const char * ); 
        void parse_quiet    ( void );
    };
}

#endif // ARGPARSE_H
