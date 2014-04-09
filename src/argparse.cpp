
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "argparse.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse
{
  const char usage[] =
  "usage: " PROGNAME " [-h] "
  "[-o OUTPUT] "
  "[-t THRESHOLD] "
  "[-a AMBIGS] "
  "[-l OVERLAP]"
  "[-d COUNTS_IN_NAME] "
  "[-f FORMAT] "
  "[-s SECOND_FASTA] "
  "[-b] "
  "[-c] "
  "[-q] "
  "[FASTA]\n";
  
  const char help_msg[] =
  "compute Tamura Nei 93 distances between aligned sequences\n"
  "\n"
  "optional arguments:\n"
  "  -h, --help               show this help message and exit\n"
  "  -o OUTPUT                direct the output to a file named OUTPUT (default=stdout)\n"
  "  -t THRESHOLD             only report (count) distances below this threshold (>=0, default=" TO_STR (DEFAULT_DISTANCE)")\n"
  "  -a AMBIGS                handle ambigous nucleotides using one of the following strategies (default=" TO_STR( DEFAULT_AMBIG ) ")\n"
  "                           resolve: resolve ambiguities to minimize distance (e.g.R matches A);\n"
  "                           average: average ambiguities (e.g.R-A is 0.5 A-A and 0.5 G-A);\n"
  "                           skip: do not include sites with ambiguous nucleotides in distance calculations;\n"
  "                           gapmm: a gap ('-') matched to anything other than another gap is like matching an N (4-fold ambig) to it;\n"
  "                           a string (e.g. RY): any ambiguity in the list is RESOLVED; any ambiguitiy NOT in the list is averaged (LIST-NOT LIST will also be averaged);\n"
  "  -f FORMAT                controls the format of the output unless -c is set (default=" TO_STR( DEFAULT_FORMAT ) ")\n"
  "                           csv: seqname1, seqname2, distance;\n"
  "                           csvn: 1, 2, distance;\n"
  "                           hyphy: {{d11,d12,..,d1n}...{dn1,dn2,...,dnn}}, where distances above THRESHOLD are set to 100;\n"
  "  -l OVERLAP               only process pairs of sequences that overlap over at least OVERLAP nucleotides (an integer >0, default=" TO_STR( DEFAULT_OVERLAP ) "):\n"
  "  -d COUNTS_IN_NAME        if sequence name is of the form 'somethingCOUNTS_IN_NAMEinteger' then treat the integer as a copy number\n"
  "                           when computing distance histograms (a character, default=" TO_STR( DEFAULT_COUNTS_IN_NAME ) "):\n"
  "  -s SECOND_FASTA          if specified, read another FASTA file from SECOND_FASTA and perform pairwise comparison BETWEEN the files (default=NULL)\n"
  "  -b                       bootstrap alignment columns before computing distances (default = false)\n"
  "                           when -s is supplied, permutes the assigment of sequences to files\n"
  "  -c                       only count the pairs below a threshold, do not write out all the pairs \n"
  "  -m                       compute inter- and intra-population means suitable for FST caclulations\n"
  "                           only applied when -s is used to provide a second file"
  "  -u PROBABILITY           subsample sequences with specified probability (a value between 0 and 1, default = " TO_STR ( DEFAULT_INCLUDE_PROB) ")\n"
  "  -q                       do not report progress updates and other diagnostics to stderr \n"
  "  FASTA                    read sequences to compare from this file (default=stdin)\n";
  
  inline
  void help()
  {
    fprintf( stderr, "%s\n%s", usage, help_msg );
    exit( 1 );
  }
  
  inline
  void ERROR( const char * msg, ... )
  {
    va_list args;
    fprintf( stderr, "%s" PROGNAME ": error: ", usage );
    va_start( args, msg );
    vfprintf( stderr, msg, args );
    va_end( args );
    fprintf( stderr, "\n" );
    exit( 1 );
  }
  
  const char * next_arg (int& i, const int argc, const char * argv[]) {
    i++;
    if (i == argc)
      ERROR ("ran out of command line arguments");
    
    return argv[i];
    
  }
  
  args_t::args_t( int argc, const char * argv[] ) :
  output( stdout ),
  input1( stdin ),
  input2( NULL ),
  distance( DEFAULT_DISTANCE ),
  ambig( DEFAULT_AMBIG ),
  format ( DEFAULT_FORMAT ),
  overlap ( DEFAULT_OVERLAP ),
  do_bootstrap( false ),
  do_count( false ),
  quiet( false ),
  do_fst( false ),
  counts_in_name ( DEFAULT_COUNTS_IN_NAME ),
  include_prob( DEFAULT_INCLUDE_PROB ),
  ambigs_to_resolve(NULL)
  {
      // skip arg[0], it's just the program name
    for (int i = 1; i < argc; ++i ) {
      const char * arg = argv[i];
      
      if ( arg[0] == '-' && arg[1] == '-' ) {
        if ( !strcmp( &arg[2], "help" ) ) help();
        else
          ERROR( "unknown argument: %s", arg );
      }
      else if ( arg[0] == '-' ) {
        if ( !strcmp( &arg[1], "h" ) ) help();
        else if (  arg[1] == 'o' ) parse_output( next_arg (i, argc, argv) );
        else if (  arg[1] == 't' ) parse_distance ( next_arg (i, argc, argv) );
        else if (  arg[1] == 'l')  parse_overlap( next_arg (i, argc, argv) );
        else if (  arg[1] == 'f')  parse_format( next_arg (i, argc, argv) );
        else if (  arg[1] == 'a')  parse_ambig( next_arg (i, argc, argv) );
        else if (  arg[1] == 's')  parse_second_in( next_arg (i, argc, argv) );
        else if (  arg[1] == 'd')  parse_counts_in_name( next_arg (i, argc, argv) );
        else if (  arg[1] == 'u')  parse_include_prob( next_arg (i, argc, argv) );
        else if (  arg[1] == 'b')  parse_bootstrap();
        else if (  arg[1] == 'c')  parse_count();
        else if (  arg[1] == 'q')  parse_quiet();
        else if (  arg[1] == 'm')  parse_fst();
        else
          ERROR( "unknown argument: %s", arg );
      }
      else
        if (i == argc-1) {
          parse_input (arg);
        } else {
          ERROR( "unknown argument: %s", arg );
        }
    }
  }
  
  args_t::~args_t() {
    if ( output && output != stdout )
      fclose( output );
    
    if ( input1 && input1 != stdin)
      fclose (input1);
    
    if ( input2 && input2 != stdin)
      fclose (input2);
  }

  
  void args_t::parse_output( const char * str )
  {
    if ( str && strcmp( str, "-" ) )
      output = fopen( str, "wb" );
    else
      output = stdout;
    
    if ( !output )
      ERROR( "failed to open the OUTPUT file %s", str );
  }
  
  void args_t::parse_input( const char * str )
  {
    if ( str && strcmp( str, "-" ) )
      input1 = fopen( str, "rb" );
    else
      input1 = stdin;
    
    if ( !input1 )
      ERROR( "failed to open the INPUT file %s", str );
    if (input1 == input2)
      ERROR( "input FASTA files must not both be stdin %s", str );
  }
  
  void args_t::parse_second_in( const char * str )
  {
    if ( str && strcmp( str, "-" ) )
      input2 = fopen( str, "rb" );
    else
      input2 = stdin;
    
    if ( !input2 )
      ERROR( "failed to open the second INPUT file %s", str );
    if (input1 == input2)
      ERROR( "input FASTA files must not both be stdin %s", str );
  }


  void args_t::parse_distance ( const char * str )
  {
    distance = atof( str );
    
    if ( distance < 0.0 || distance > 1.0)
      ERROR( "genetic distance threshold must be in [0,1], had: %s", str );
  }

  void args_t::parse_include_prob ( const char * str )
  {
    include_prob = atof( str );
    
    if ( include_prob < 0.0 || include_prob > 1.0)
      ERROR( "sequence inclusion probability must be in [0,1], had: %s", str );
  }

  void args_t::parse_counts_in_name ( const char * str )
  {
    counts_in_name = str[0];
    
    if ( ! isprint (counts_in_name))
      ERROR( "count separator must be a printable character, had: %s", str );
  }

  void args_t::parse_overlap ( const char * str )
  {
    overlap = atoi( str );
    
    if ( overlap == 0 )
      ERROR( "overlap must be positive, had: %s", str );
  }
  
  void args_t::parse_ambig( const char * str )
  {
    if (!strcmp (str, "resolve")) {
      ambig = resolve;
    } else if (!strcmp (str, "average")) {
      ambig = average;
    } else if (!strcmp (str, "skip")) {
      ambig = skip;
    } else if (!strcmp (str, "gapmm")) {
      ambig = gapmm;
    } else {
      ambig = subset;
      ambigs_to_resolve = new char [strlen (str) + 1];
      strcpy (ambigs_to_resolve, str);
    }
  }
  
  void args_t::parse_format( const char * str )
  {
    if (!strcmp (str, "csv")) {
      format = csv;
    } else if (!strcmp (str, "csvn")) {
      format = csvn;
    } else if (!strcmp (str, "hyphy")) {
      format = hyphy;
    } else  {
      ERROR( "invalid output format: %s", str );
    }
  }
  
  void args_t::parse_count()
  {
    do_count = true;
  }
  
  void args_t::parse_bootstrap()
  {
    do_bootstrap = true;
  }
  
  void args_t::parse_quiet()
  {
    quiet = true;
  }

  void args_t::parse_fst()
  {
    do_fst = true;
  }
}
