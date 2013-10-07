
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "argparse_trim.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse
{
  const char usage[] =
  "usage: " PROGNAME " [-h] "
  "[-o OUTPUT] "
  "[-a AMBIGS] "
  "[-s START] "
  "[-e END] "
  "[-c COVERAGE] "
  "[-t DATATYPE] "
  "[-q] "
  "[FASTA]\n";
  
  const char help_msg[] =
  "read a FASTA MSA and output only the sequences which overlap a given coordinate range\n"
  "\n"
  "optional arguments:\n"
  "  -h, --help               show this help message and exit\n"
  "  -o OUTPUT                direct the FASTA file with matching (and trimmed reads) to a file named OUTPUT (default=stdout)\n"
  "  -a AMBIGS                count the following characters AGAINST coverage numbers (default=" TO_STR( DEFAULT_AMBIG ) ")\n"
  "                           gaponly: gap character ('-') is counted as lack of coverage;\n"
  "                           nfold: gap character AND N-fold ambiguities (N and ?) are counted as lack of coverage;\n"
  "                           threefold: gap character, 4-fold ambiguities, AND 3-fold ambiguities (e.g. M and S) are counted as lack of coverage;\n"
  "                           any: ALL incompletely resolved characters are counted as lack of coverage;\n"
  "  -s START                 start of the region to filter, 0-based, INCLUSIVE, must be an integer strictly than the length of the alignment (default=" TO_STR( DEFAULT_START ) ")\n"
  "                           e.g. -s 102 -e 203 will span nucleotides 103 through 204\n"
  "  -e END                   end of the region to filter, INCLUSIVE, 0-based, capped at [length of the alignment - 1] (default=" TO_STR( DEFAULT_END ) ")\n"
  "  -c COVERAGE              require that retained reads cover at least this proportion of the region (default=" TO_STR( DEFAULT_COVERAGE ) ")\n"
  "                           must be a floating point in (0,1]\n"
  "  -t DATATYPE              the type of data expected (default=" TO_STR( DEFAULT_DATA ) ")\n"
  "                           dna: DNA or RNA (IUPAC);\n"
  "                           protein : protein (IUPAC);\n"
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
