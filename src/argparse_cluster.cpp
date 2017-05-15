
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "argparse_cluster.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse
{
  const char usage[] =
  "usage: " PROGNAME " [-h] "
  "[-o OUTPUT] "
  "[-a AMBIGS] "
  "[-c COVERAGE] "
  "[-t DATATYPE] "
  "[-q] "
  "[FASTA]\n";
  
  const char help_msg[] =
  "read a FASTA MSA and (greedy) cluster sequences that lie within a specific distance of each other\n"
  "\n"
  "optional arguments:\n"
  "  -h, --help               show this help message and exit\n"
  "  -o OUTPUT                direct the output file with clusters to OUTPUT either \n"
  "                           JSON or sets of FASTA files representing individual clusters\n"
  "                           (default=stdout) see also -f\n"
  "  -t THRESHOLD             sequences which lie within this distance will be clustered (>=0, default=" TO_STR (DEFAULT_DISTANCE)")\n"
  "  -a AMBIGS                handle ambigous nucleotides using one of the following strategies (default=" TO_STR( DEFAULT_AMBIG ) ")\n"
  "                           resolve: resolve ambiguities to minimize distance (e.g.R matches A);\n"
  "                           average: average ambiguities (e.g.R-A is 0.5 A-A and 0.5 G-A);\n"
  "                           skip: do not include sites with ambiguous nucleotides in distance calculations;\n"
  "                           gapmm: a gap ('-') matched to anything other than another gap is like matching an N (4-fold ambig) to it;\n"
  "                           a string (e.g. RY): any ambiguity in the list is RESOLVED; any ambiguitiy NOT in the list is averaged \n"
  "                           (LIST-NOT LIST will also be averaged);\n"
  "  -c CLUSTER-TYPE          create clusters based on the following rules\n"
  "     "
  
  
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
  start_coord ( DEFAULT_START ),
  end_coord ( DEFAULT_END ),
  coverage ( DEFAULT_COVERAGE ),
  data ( DEFAULT_DATA ), 
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
        else if (  arg[1] == 's')  parse_start ( next_arg (i, argc, argv) );
        else if (  arg[1] == 'e')  parse_end( next_arg (i, argc, argv) );
        else if (  arg[1] == 'c')  parse_coverage ( next_arg (i, argc, argv) );
        else if (  arg[1] == 't')  parse_data ( next_arg (i, argc, argv) );
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
  

  void args_t::parse_start ( const char * str )
  {
    start_coord = atoi (str);
    
    if ( start_coord > end_coord )
      ERROR( "start of the filtering frame must be less than end of the frame, had: %ld - %ld", start_coord, end_coord );
  }

  void args_t::parse_end ( const char * str )
  {
    end_coord = atoi (str);
    
    if ( start_coord > end_coord )
      ERROR( "start of the filtering frame must be less than end of the frame, had: %ld - %ld", start_coord, end_coord );
  }
  
  void args_t::parse_coverage ( const char * str )
  {
    coverage = atof (str);
    
    if ( coverage <= 0.0 || coverage > 1.0 )
      ERROR( "coverage must be in (0,1], had: %g", coverage );
  }
  
      void args_t::parse_ambig( const char * str )
  {
    if (!strcmp (str, "gaponly")) {
      ambig = gaponly;
    } else if (!strcmp (str, "fourfold")) {
      ambig = nfold;
    } else if (!strcmp (str, "threefold")) {
      ambig = threefold;
    } else if (!strcmp (str, "any")) {
      ambig = any;
    } else  {
      ERROR( "invalid ambiguity mode type: %s", str );
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
