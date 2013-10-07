
#ifndef ARGPARSE_H
#define ARGPARSE_H

// argument defaults


#define PROGNAME                 "readreduce"
#define DEFAULT_AMBIG             resolve
#define DEFAULT_COUNTS_IN_NAME    ':'
#define DEFAULT_OVERLAP           100
#define DEFAULT_SIZE              16


namespace argparse
{
  

    enum ambig_t {
      resolve,
      average,
      skip,
      gapmm
    };

  class args_t
    {
    public:
 
        FILE            * output,
                        * input;
                                     
        ambig_t         ambig;
        unsigned long   overlap,
                        cluster_size;
                        
        bool            quiet,
                        json;
                        
        char            counts_in_name;
        
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input    ( const char * );
        void parse_output    ( const char * );
        void parse_overlap  ( const char * );
        void parse_size     ( const char * );
        void parse_ambig    ( const char * );
        void parse_counts_in_name   ( const char * );
        void parse_quiet ( void );
        void parse_json ( void );
      
    };
}

#endif // ARGPARSE_H
