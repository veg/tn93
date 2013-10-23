
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "argparse_cf.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse
{
  const char usage[] =
  "usage: " PROGNAME " [-h] "
  "[-o OUTPUT] "
  "[-a AMBIGS] "
  "[-t TYPE] "
  "[-d COUNTS_IN_NAME] "
  "[-q] "
  "[FASTA]\n";
  
  const char help_msg[] =
  "compute Tamura Nei 93 distances between aligned sequences\n"
  "\n"
  "optional arguments:\n"
  "  -h, --help               show this help message and exit\n"
  "  -o OUTPUT                direct the output JSON to a file named OUTPUT (default=stdout)\n"
  "  -a AMBIGS                handle ambigous characters using one of the following strategies (default=" TO_STR( DEFAULT_AMBIG ) ")\n"
  "                           average: average ambiguities (e.g. a nucleotide R adds 0.5 to A and G coverage for that position);\n"
  "                           ignore: ambiguities contribute nothing to coverage counts;\n"
  "  -t DATATYPE              the type of data expected (default=" TO_STR( DEFAULT_FORMAT ) ")\n"
  "                           dna: DNA or RNA (IUPAC);\n"
  "                           protein : protein (IUPAC);\n"
  "  -d COUNTS_IN_NAME        if sequence name is of the form 'somethingCOUNTS_IN_NAMEinteger' then treat the integer as a copy number\n"
  "                           when computing coverages (a character, default=" TO_STR( COUNTS_IN_NAME ) "):\n"
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
  data ( DEFAULT_DATA ),
  counts_in_name ( DEFAULT_COUNTS_IN_NAME ),
  quiet (false)
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
        else if (  arg[1] == 'a')  parse_ambig( next_arg (i, argc, argv) );
        else if (  arg[1] == 't')  parse_data ( next_arg (i, argc, argv) );
        else if (  arg[1] == 'd')  parse_counts_in_name( next_arg (i, argc, argv) );
        else if (  arg[1] == 'q')  parse_quiet ( );
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
    
    if ( output == NULL )
      ERROR( "failed to open the OUTPUT file %s", str );
  }
  
  void args_t::parse_input( const char * str )
  {
    if ( str && strcmp( str, "-" ) )
      input = fopen( str, "rb" );
    else
      input = stdin;
    
    if ( input == NULL )
      ERROR( "failed to open the INPUT file %s", str );
  }
  

  void args_t::parse_counts_in_name ( const char * str )
  {
    counts_in_name = str[0];
    
    if ( ! isprint (counts_in_name))
      ERROR( "count separator must be a printable character, had: %s", str );
  }

  void args_t::parse_ambig( const char * str )
  {
    if (!strcmp (str, "average")) {
      ambig = average;
    } else if (!strcmp (str, "ignore")) {
      ambig = ignore;
    } else {
      ERROR( "invalid ambiguity handling mode: %s", str );
    }
  }
  
  void args_t::parse_data( const char * str )
  {
    if (!strcmp (str, "dna")) {
      data = dna;
    } else if (!strcmp (str, "protein")) {
      data = protein;
    } else  {
      ERROR( "invalid data type: %s", str );
    }
  }
  
  void args_t::parse_quiet ( void ) {
    quiet = true;
  }
}
