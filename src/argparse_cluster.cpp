
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
    "[-l OVERLAP] "
    "[-t THERSHOLD] "
    "[-c CLUSTER-TYPE] "
    "[-m OUTPUT_MODE] "
    "[-g FRACTION] "
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
    "                           skip:    do not include sites with ambiguous nucleotides in distance calculations;\n"
    "                           gapmm:   a gap ('-') matched to anything other than another gap is like matching an N (4-fold ambig) to it;\n"
    "                           a string (e.g. RY): any ambiguity in the list is RESOLVED; any ambiguitiy NOT in the list is averaged \n"
    "                           (LIST-NOT LIST will also be averaged);\n"
    "  -c CLUSTER-TYPE          create clusters based on the following rules (default=" TO_STR (DEFAULT_CLUSTER_TYPE) ")\n"
    "                           all:     each sequence in the cluster is within the specified distance threshold \n"
    "                                    of EVERY other sequence; this cluster definition is stable from run to run\n"
    "                           any:     each sequence in the cluster is within the specified distance threshold \n"
    "                                    of AT LEAST ONE other sequence; this cluster definition need NOT be stable \n"
    "                                    from run to run\n"
    "  -m OUTPUT-MODE           output clusters in the following format (default=" TO_STR (DEFAULT_OUTPUT_MODE) ")\n"
    "                           json:    write a JSON file with individual clusters represented by the 'centroid' \n"
    "                                    i.e. the longest sequence, and the list of all other sequence members\n"
    "                           files:   if OUTPUT is set then write each cluster to as a FASTA file OUTPUT.x where\n"
    "                                    x is the cluster ID; if OUTPUT is stdout, then write each cluster as a FASTA\n"
    "                                    file separated by a line of ======\n"
    "  -l OVERLAP               only process pairs of sequences that overlap over at least OVERLAP nucleotides \n"
    "                           (an integer >0, default=" TO_STR( DEFAULT_OVERLAP ) ")\n"
    "  -g FRACTION              in combination with AMBIGS, works to limit (for resolve and string options to AMBIG)\n"
    "                           the maximum tolerated FRACTION of ambiguous characters; sequences whose pairwise comparisons\n"
    "                           include no more than FRACTION [0,1] of sites with resolvable ambiguities will be resolved\n"
    "                           while all others will be AVERAGED (default = " TO_STR ( DEFAULT_FRACTION ) ")\n"
    "  -f FIRST                 treat first sequence as regular sequence (default is to treat it as reference and skip)\n"
    "  -q QUIET                 do not print progress updates to stderr (default is to print)\n";

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
    ambig( DEFAULT_AMBIG ),
    cluster_type ( DEFAULT_CLUSTER_TYPE ),
    output_mode ( DEFAULT_OUTPUT_MODE ),
    overlap ( DEFAULT_OVERLAP ),
    quiet (false),
    ambigs_to_resolve (NULL),
    trunk_path (NULL),
    resolve_fraction ( DEFAULT_FRACTION ),
    first_regular (false) {
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
                else if (  arg[1] == 't')  parse_distance ( next_arg (i, argc, argv) );
                else if (  arg[1] == 'c')  parse_cluster ( next_arg (i, argc, argv) );
                else if (  arg[1] == 'm')  parse_output_mode ( next_arg (i, argc, argv) );
                else if (  arg[1] == 'g')  parse_fraction ( next_arg (i, argc, argv) );
                else if (  arg[1] == 'l')  parse_overlap ( next_arg (i, argc, argv) );
                else if (  arg[1] == 'q')  parse_quiet ( );
                else if (  arg[1] == 'f')  parse_first ( );
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

        if (ambigs_to_resolve){
            delete [] ambigs_to_resolve;
        }

        if (trunk_path) {
            delete [] trunk_path;
        }

        if ( input && input != stdin)
            fclose (input);

    }


    void args_t::parse_output( const char * str ) {
        if ( str && strcmp( str, "-" ) ) {
            trunk_path = new char [strlen (str) + 1];
            strcpy (trunk_path, str);
        }
    }

    void args_t::parse_distance ( const char * str )
    {
        distance = atof( str );

        if ( distance < 0.0 || distance > 1.0)
            ERROR( "genetic distance threshold must be in [0,1], had: %s", str );
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

    void args_t::parse_cluster( const char * str ) {
        if (!strcmp (str, "all")) {
            cluster_type = all;
        } else if (!strcmp (str, "any")) {
            cluster_type = any;
        } else  {
            ERROR( "invalid cluster construction mode type: %s", str );
        }
    }

    void args_t::parse_output_mode( const char * str ) {
        if (!strcmp (str, "json")) {
            output_mode = json;
        } else if (!strcmp (str, "files")) {
            output_mode = files;
        } else  {
            ERROR( "invalid output format mode type: %s", str );
        }
    }

    void args_t::parse_fraction ( const char * str ) {
        resolve_fraction = atof( str );
        if ( resolve_fraction < 0.0 || resolve_fraction > 1.0)
            ERROR( "resolve ambigous fraction must be in [0,1], had: %s", str );
    }

    void args_t::parse_overlap ( const char * str ){
        overlap = atoi( str );

        if ( overlap == 0L )
            ERROR( "overlap must be positive, had: %s", str );
    }

    void args_t::parse_ambig( const char * str ) {
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


    void args_t::parse_quiet ( void ) {
        quiet = true;
    }

    void args_t::parse_first ( void ) {
        first_regular = true;
    }

}
