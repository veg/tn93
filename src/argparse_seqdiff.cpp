
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "argparse_seqdiff.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse
{
  const char usage[] =
  "usage: " PROGNAME " [-h] "
  "[-o OUTPUT] "
  "[-a AMBIGS] "
  "[-f FORMAT] "
  "[-q] "
  "[-c] "
  "-r REFERENCE "
  "[FASTA]\n";
  
  const char help_msg[] =
  "read a FASTA MSA and a reference sequence and output \n"
  "\n"
  "optional arguments:\n"
  "  -h, --help               show this help message and exit\n"
  "  -o OUTPUT                direct the file with recorded differences from reference to a file named OUTPUT (default=stdout)\n"
  "  -r REFERENCE             requried file for the reference sequence (i.e. the sequence from which the differences are computed)\n"
  "  -f FORMAT                controls the format of the output (default=" TO_STR( DEFAULT_FORMAT ) ")\n"
  "                           json: compressed JSON format; 'variants' : [list of variants], 'sequences' : [list of variant codes]\n"
  "                           csv: sequence_id, position, alternative nucleotide;\n"
  "  -a AMBIGS                when comparing to reference tree ambiguities as follows (default=" TO_STR( DEFAULT_AMBIG ) ")\n"
  "                           informative: compare characters (i.e. different character => mismatch), except for N -- a 4-fold misisng base;\n"
  "                           resolve: if an ambiguity is resolvable to match the reference, do not count it;\n"
  "                           all: report all differences;\n"
  "  -q                       do not report progress updates and other diagnostics to stderr \n"
  "  -c                       DO NOT combine all identical sequences into a group (default is TO COMBINE); only applies to the JSON format\n"
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
  reference (NULL),
  ambig( DEFAULT_AMBIG ),
  out_format ( DEFULT_FORMAT ),
  quiet (false),
  compress (true)
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
        else if (  arg[1] == 'r')  parse_reference ( next_arg (i, argc, argv) );
        else if (  arg[1] == 'f')  parse_format( next_arg (i, argc, argv) );
        else if (  arg[1] == 'q')  parse_quiet ( );
        else if (  arg[1] == 'c')  parse_compress ( );
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
      
      if (!reference) {
          ERROR( "missing required argument: %s", "-r" );
      }
  }
  
  args_t::~args_t() {
    if ( output && output != stdout )
      fclose( output );
    
    if ( input && input != stdin)
      fclose (input);
      
    if ( reference && reference != stdin)
     fclose (reference);
    
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
  

  void args_t::parse_reference( const char * str )
  {
      reference = fopen( str, "rb" );
      
      if ( reference == NULL )
        ERROR( "failed to open the REFERENCE file %s", str );
  }

  
  
  void args_t::parse_ambig( const char * str )
  {
    if (!strcmp (str, "all")) {
      ambig = all;
    } else if (!strcmp (str, "resolve")) {
      ambig = resolve;
    } else if (!strcmp (str, "informative")) {
        ambig = informative;
    } else  {
      ERROR( "invalid ambiguity mode type: %s", str );
    }
  }
  
  void args_t::parse_format( const char * str )
  {
    if (!strcmp (str, "json")) {
       out_format = json;
    } else if (!strcmp (str, "csv")) {
       out_format = csv;
    } else  {
      ERROR( "invalid output format: %s", str );
    }
  }

  void args_t::parse_quiet ( void ) {
    quiet = true;
  }

  void args_t::parse_compress ( void ) {
    compress = false;
  }
}
