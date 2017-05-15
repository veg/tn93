
#ifndef ARGPARSE_H
#define ARGPARSE_H

// argument defaults

#define PROGNAME                 "tn93-cluster"
#define DEFAULT_AMBIG             resolve
#define DEFAULT_FRACTION          1.0
#define DEFAULT_DISTANCE          0.015
#define DEFAULT_OVERLAP           100
#define DEFAULT_CLUSTER_OP        all
#define DEFAULT_OUTPUT_MODE       json

#ifndef VERSION_NUMBER
#define VERSION_NUMBER            "UNKNOWN"
#endif

namespace argparse {
  
  enum ambig_t {
    resolve,
    average,
    skip,
    gapmm,
    subset
  };
  
  enum cluster_t {
    all,
    any
  };
  
  enum output_t {
    json,
    files
  };
  
  class args_t {
  public:
    
    FILE            * output,
                    * input;
    
    double          distance;
    ambig_t         ambig;
    cluster_t       cluster_mode;
    unsigned long   overlap;
    bool            quiet;
    char            *ambigs_to_resolve;
    double          resolve_fraction;
    
    args_t( int, const char ** );
    ~args_t();
    
  private:
    void parse_input    ( const char * );
    void parse_output   ( const char * );
    void parse_distance ( const char * );
    void parse_overlap  ( const char * );
    void parse_ambig    ( const char * );
    void parse_quiet    ( void );
    void parse_fraction ( const char *);
    void parse_cluster  ( const char *);
    
  };
}

#endif // ARGPARSE_H
