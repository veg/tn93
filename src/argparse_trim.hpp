
#ifndef ARGPARSE_H
#define ARGPARSE_H
#include <cstdarg>

// argument defaults


#define PROGNAME                 "selectreads"

#define        DEFAULT_AMBIG         gaponly
#define        DEFAULT_START         0L
#define        DEFAULT_END           0xFFFFFF
#define        DEFAULT_COVERAGE      0.95
#define        DEFAULT_DATA          dna


namespace argparse
{
  void ERROR( const char * msg, ... );


  enum ambig_mode {
      gaponly,
      nfold,
      threefold,
      any
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
             
        ambig_mode      ambig;
        unsigned        long start_coord;
        unsigned        long end_coord;
        double          coverage;
        data_t          data;
        bool            quiet;
       
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input    ( const char * );
        void parse_output   ( const char * );
        void parse_quiet    ( void );
        void parse_start    ( const char * );
        void parse_end      ( const char * );
        void parse_coverage ( const char * );
        void parse_data     ( const char * );  
        void parse_ambig    ( const char * );
    };
}

#endif // ARGPARSE_H
