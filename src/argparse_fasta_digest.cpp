
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "argparse_fasta_digest.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse {
    const char usage[] =
    "usage: " PROGNAME " [-h] "
    "[-o OUTPUT] "
    "[-r REVERSE_COMPLEMENT] "
    "[-m MATCH] "
    "-e MOTIFS "
    "[FASTA]\n";
    
    const char help_msg[] =
    "Read a FASTA file, a comma separated list of nucleotide binding motifs and report, for each sequence, 1-based coordinates of any motifs (nucleotide strings, 4-16 letters) found in the sequence\n"
    "\n"
    "optional arguments:\n"
    "  -h, --help               show this help message and exit\n"
    "  -o OUTPUT                write the resulting JSON file to a file named OUTPUT (default=stdout)\n"
    "  -r OPERATION             consider reverse complements (default=" TO_STR( DEFAULT_RC ) ")\n"
    "                           complement: for each motif, also consider its reverse complement  \n"
    "                           nocomplement: only consider the motifs themselves\n"
    "  -m MATCH                 how to match the motif (default=" TO_STR( DEFAULT_MATCH ) ")\n"
    "                           exact: match the motif exactly\n"
    "                           oneoff: allow any one motif mis-match\n"
    "                           transition: allow any one motif mis-match that is a transition (A:G, C:T)\n"
    "  -e MOTIF_LIST            specify delimited nucleotide sequences to use as motifs\n"
    "                           not case sensitve (aaaccc is the same as AAACCC)\n"
    "                           any contigous set of ACGTs define a motif but it must be 4-12 bases in length\n"
    "  FASTA                    read FASTA sequences from this file (default=stdin)\n";
    
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
    input( stdin ),
    output( stdout ),
    rc ( DEFAULT_RC ),
    match (DEFAULT_MATCH),
    motif_list (NULL)
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
                else if (  arg[1] == 'r')  parse_rc_option ( next_arg (i, argc, argv) );
                else if (  arg[1] == 'm')  parse_match_mode( next_arg (i, argc, argv) );
                else if (  arg[1] == 'e')  parse_motif( next_arg (i, argc, argv) );
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
                
        if (motif_list == NULL) {
            ERROR ("Required argument MOTIF_LIST was not provided");
        }
    }
    
    args_t::~args_t() {
        if ( output && output != stdout )
            fclose( output );
        
        if ( input && input != stdin)
            fclose ( input );
        
        if ( motif_list )
            free ( motif_list );
    }
    
    
    void args_t::parse_output( const char * str ) {
        if ( str && strcmp( str, "-" ) )
            output = fopen( str, "wb" );
        else
            output = stdout;
        
        if ( output == NULL )
            ERROR( "failed to open the OUTPUT file %s", str );
    }
    
    void args_t::parse_input ( const char * str ) {
        if ( str && strcmp( str, "-" ) )
            input = fopen( str, "rb" );
        else
            input = stdin;
        
        if ( input == NULL )
            ERROR( "failed to open the INPUT file %s", str );
    }
    
    
    
    void args_t::parse_rc_option( const char * str )
    {
        if (!strcmp (str, "complement")) {
            rc = complement;
        } else if (!strcmp (str, "nocomplement")) {
            rc = nocomplement;
        } else  {
            ERROR( "invalid file operation: %s", str );
        }
    }
    
    void args_t::parse_match_mode ( const char * str )
    {
        if (!strcmp (str, "exact")) {
            match = exact;
        } else if (!strcmp (str, "oneoff")) {
            match = oneoff;
        } else if (!strcmp (str, "transitions")) {
            match = transition;
        } else  {
            ERROR( "invalid match mode: %s", str );
        }
    }
    
    void args_t::parse_motif ( const char * str  ) {
        const int L = strlen (str);
        if (L >= 4) {
            motif_list = (char*)malloc (L + 1);
            strcpy (motif_list, str);
        } else {
            ERROR( "motif list must be at least 4 characters long: %s", str );
        }
    }
}
