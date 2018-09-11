
#ifndef ARGPARSE_H
#define ARGPARSE_H

// argument defaults

#define PROGNAME                 "tn93-cluster"
#define DEFAULT_AMBIG             resolve
#define DEFAULT_FRACTION          1.0
#define DEFAULT_DISTANCE          0.015
#define DEFAULT_OVERLAP           100
#define DEFAULT_CLUSTER_TYPE      all
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

    FILE            * input;

    double          distance;
    ambig_t         ambig;
    cluster_t       cluster_type;
    output_t        output_mode;
    unsigned long   overlap;
    bool            quiet;
    char            *ambigs_to_resolve;
    char            *trunk_path;
    double          resolve_fraction;
    bool            first_regular;

    args_t( int, const char ** );
    ~args_t();

  private:
    void parse_input    ( const char * );
    void parse_output   ( const char * );
    void parse_output_mode   ( const char * );
    void parse_distance ( const char * );
    void parse_overlap  ( const char * );
    void parse_ambig    ( const char * );
    void parse_quiet    ( void );
    void parse_fraction ( const char *);
    void parse_cluster  ( const char *);
    void parse_first    (void);

  };
}

#endif // ARGPARSE_H
