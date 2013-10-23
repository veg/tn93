
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "argparse_merge.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse
{
  const char usage[] =
  "usage: " PROGNAME " [-h] "
  "[-o OUTPUT] "
  "[-a AMBIGS] "
  "[-l OVERLAP]"
  "[-d COUNTS_IN_NAME] "
  "[-s SIZE] "
  "[-q] "
  "[-j] "
  "[FASTA]\n";
  
  const char help_msg[] =
  "merge matching (subject to ambig resolution) reads into clusters\n"
  "\n"
  "optional arguments:\n"
  "  -h, --help               show this help message and exit\n"
  "  -o OUTPUT                direct the output to a file named OUTPUT (default=stdout)\n"
  "  -a AMBIGS                handle ambigous nucleotides using one of the following strategies (default=" TO_STR( DEFAULT_AMBIG ) ")\n"
  "                           resolve: resolve ambiguities to minimize distance (e.g.R matches A);\n"
  "                           average: average ambiguities (e.g.R-A is 0.5 A-A and 0.5 G-A);\n"
  "                           skip: do not include sites with ambiguous nucleotides in distance calculations;\n"
  "                           gapmm: a gap ('-') matched to anything other than another gap is like matching an N (4-fold ambig) to it;\n"
  "  -l OVERLAP               merge reads that match over at least this many positions (an integer >0, default=" TO_STR( DEFAULT_OVERLAP ) "):\n"
  "  -s SIZE                  merge reads that match over at least this many positions (an integer >0, default=" TO_STR( DEFAULT_SIZE ) "):\n"
  "  -d COUNTS_IN_NAME        if sequence name is of the form 'somethingCOUNTS_IN_NAMEinteger' then treat the integer as a copy number\n"
  "                           when counting reads; also output cluster sizes using the same separator (a character, default=" TO_STR( DEFAULT_COUNTS_IN_NAME ) "):\n"
  "  -j                       report the results to a JSON file (instead of a FASTA MSA) \n"
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
  input( stdin ),
  ambig( DEFAULT_AMBIG ),
  overlap ( DEFAULT_OVERLAP ),
  cluster_size ( DEFAULT_SIZE ),
  quiet( false ),
  json ( false ),
  counts_in_name ( DEFAULT_COUNTS_IN_NAME )
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
        else if (  arg[1] == 'l')  parse_overlap( next_arg (i, argc, argv) );
        else if (  arg[1] == 's')  parse_size( next_arg (i, argc, argv) );
        else if (  arg[1] == 'a')  parse_ambig( next_arg (i, argc, argv) );
        else if (  arg[1] == 'd')  parse_counts_in_name ( next_arg (i, argc, argv) );
        else if (  arg[1] == 'q')  parse_quiet();
        else if (  arg[1] == 'j')  parse_json();
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
    
    if ( input && input != stdin)
      fclose (input);
    
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
      input = fopen( str, "rb" );
    else
      input = stdin;
    
    if ( !input )
      ERROR( "failed to open the INPUT file %s", str );
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
    
    if ( overlap == 0L )
      ERROR( "overlap must be positive, had: %s", str );
  }
  
  void args_t::parse_size ( const char * str )
  {
    cluster_size = atoi( str );
    
    if ( cluster_size == 0L )
      ERROR( "cluster size must be positive, had: %s", str );
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
      ERROR( "invalid ambiguity handling mode: %s", str );
    }
  }
  
  void args_t::parse_quiet()
  {
    quiet = true;
  } 
  
  void args_t::parse_json()
  {
    json = true;
  }    
}
