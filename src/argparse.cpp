
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "argparse.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse
{
    const char usage[] =
        "usage: " tn93 " [-h] "
        "[-o OUTPUT] "
        "[-t THRESHOLD] "
        "[-a AMBIGS] "
        "[-l OVERLAP]"
        "[-d COUNTS_IN_NAME] "
        "[-f FORMAT] "
        "[-s SECOND_FASTA] "
        "[-b] "
        "[-c] "
        "[FASTA]\n";

    const char help_msg[] =
        "compute Tamura Nei 93 distances between aligned sequences\n"
        "\n"
        "optional arguments:\n"
        "  -h, --help               show this help message and exit\n"
        "  -o OUTPUT                direct the output to a file named OUTPUT (default=stdout)\n"
        "  -t THRESHOLD             only report (count) distances below this threshold (>=0, default=" TO_STR (DEFAULT_DISTANCE)")\n"
        "  -a AMBIGS                handle ambigous nucleotides using one of the following strategies (default=" TO_STR( DEFAULT_AMBIG ) ")\n"
        "                           RESOLVE: resolve ambiguities to minimize distance (e.g.R matches A);\n"
        "                           AVERAGE: average ambiguities (e.g.R-A is 0.5 A-A and 0.5 G-A);\n"
        "                           SKIP: do not include sites with ambiguous nucleotides in distance calculations;\n"
        "                           GAPMM: a gap ('-') matched to anything other than another gap is like matching an N (4-fold ambig) to it;\n"
        "  -l OVERLAP               only process pairs of sequences that overlap over at least OVERLAP nucleotides (an integer >=0, default=" TO_STR( DEFAULT_OVERLAP ) "):\n"
        "  -d COUNTS_IN_NAME        if sequence name is of the form 'somethingCOUNTS_IN_NAMEinteger' then treat the integer as a copy number\n"
        "                           when computing distance histograms (a character, default=" TO_STR( COUNTS_IN_NAME ) "):\n"
        "  -s SECOND_FASTA          if specified, read another FASTA file from SECOND_FASTA and perform pairwise comparison BETWEEN the files (default=NULL)\n"
        "  -b                       bootstrap alignment columns before computing distances (default = false)\n"
        "  -c                       only count the pairs below a threshold, do not write out all the pairs \n"
        "  FASTA                    read sequences to compare from this file (default=stdin)n";

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
        fprintf( stderr, "%s" QFILT ": error: ", usage );
        va_start( args, msg );
        vfprintf( stderr, msg, args );
        va_end( args );
        fprintf( stderr, "\n" );
        exit( 1 );
    }

    args_t::args_t( int argc, const char * argv[] ) :
        fasta( NULL ),
        fastq( NULL ),
        qual( NULL ),
        output( stdout ),
        min_length( DEFAULT_MIN_LENGTH ),
        min_qscore( DEFAULT_MIN_QSCORE ),
        json( false ),
        punch( '\0' ),
        tag_length( 0 ),
        tag_mismatch( DEFAULT_TAG_MISMATCH ),
        format( DEFAULT_FORMAT )
    {
        int i;
        // make sure tag is an empty string
        tag[0] = '\0';
        // handle the mode separately
        parse_mode( TO_STR( DEFAULT_MODE ) );

        // skip arg[0], it's just the program name
        for ( i = 1; i < argc; ++i ) {
            const char * arg = argv[i];

            if ( arg[0] == '-' && arg[1] == '-' ) {
                if ( !strcmp( &arg[2], "help" ) ) help();
#if 0
                else if ( !strcmp( &arg[2], "fastq" ) ) parse_fastq( argv[++i] );
                else if ( !strcmp( &arg[2], "fasta" ) ) {
                    parse_fasta( argv[i + 1], argv[i + 2] );
                    i += 2;
                }
                else if ( !strcmp( &arg[2], "minlength" ) ) parse_minlength( argv[++i] );
                else if ( !strcmp( &arg[2], "minqscore" ) ) parse_minqscore( argv[++i] );
                else if ( !strcmp( &arg[2], "mode" ) ) parse_mode( argv[++i] );
                else if ( !strcmp( &arg[2], "tag" ) ) parse_tag( argv[++i] );
                else if ( !strcmp( &arg[2], "tagmismatch" ) ) parse_tagmismatch( argv[++i] );
#endif
                else
                    ERROR( "unknown argument: %s", arg );
            }
            else if ( arg[0] == '-' ) {
                if ( !strcmp( &arg[1], "h" ) ) help();
                else if ( !strcmp( &arg[1], "F" ) ) {
                    parse_fasta( argv[i + 1], argv[i + 2] );
                    i += 2;
                }
                else if ( !strcmp( &arg[1], "Q" ) ) parse_fastq( argv[++i] );
                else if ( !strcmp( &arg[1], "o" ) ) parse_output( argv[++i] );
                else if ( !strcmp( &arg[1], "l" ) ) parse_minlength( argv[++i] );
                else if ( !strcmp( &arg[1], "q" ) ) parse_minqscore( argv[++i] );
                else if ( !strcmp( &arg[1], "m" ) ) parse_mode( argv[++i] );
                else if ( !strcmp( &arg[1], "s" ) ) parse_split();
                else if ( !strcmp( &arg[1], "p" ) ) parse_hpoly();
                else if ( !strcmp( &arg[1], "a" ) ) parse_ambig();
                else if ( !strcmp( &arg[1], "j" ) ) parse_json();
                else if ( !strcmp( &arg[1], "P" ) ) parse_punch( argv[++i] );
                else if ( !strcmp( &arg[1], "T" ) ) parse_tag( argv[++i] );
                else if ( !strcmp( &arg[1], "t" ) ) parse_tagmismatch( argv[++i] );
                else if ( !strcmp( &arg[1], "f" ) ) parse_format( argv[++i] );
                else
                    ERROR( "unknown argument: %s", arg );
            }
            else
                ERROR( "unknown argument: %s", arg );
        }

        if ( !fastq && ( !fasta || !qual ) )
            ERROR( "missing required argument -F FASTA QUAL or -Q FASTQ" );

        if ( punch && ( split || hpoly || ambig ) )
            ERROR( "-P CHAR is incompatible with any of -s, -p, and -a" );
    }

    args_t::~args_t() {
        if ( fasta )
            delete fasta;
        if ( fastq )
            delete fastq;
        if ( qual )
            delete qual;
        if ( output && output != stdin )
            fclose( output );
    }

    void args_t::parse_fasta( const char * fstr, const char * qstr )
    {
        if ( fastq )
            ERROR( "-F and -Q are mutually exclusive" );

        if ( !strcmp( fstr, "-" ) && !strcmp( qstr, "-" ) )
            ERROR( "only one argument to -F FASTA and QUAL can be STDIN" );

        fasta = new ifile::ifile_t( fstr );
        qual = new ifile::ifile_t( qstr );

        if ( !fasta || !fasta->good() )
            ERROR( "failed to open the FASTA file %s", fstr );

        if ( !qual && !qual->good() )
            ERROR( "failed to open the QUAL file %s", qstr );
    }

    void args_t::parse_fastq( const char * str )
    {
        if ( fasta || qual )
            ERROR( "-Q and -F are mutually exclusive" );

        fastq = new ifile::ifile_t( str );

        if ( !fastq->good() )
            ERROR( "failed to open the FASTQ file %s", str );
    }

    void args_t::parse_output( const char * str )
    {
        if ( str && strcmp( str, "-" ) )
            output = fopen( str, "wb" );
        else
            output = stdin;

        if ( !output )
            ERROR( "failed to open the OUTPUT file %s", str );
    }

    void args_t::parse_minlength( const char * str )
    {
        long val = atoi( str );

        if ( val < 1 )
            ERROR( "minimum length expected a positive integer, had: %s", str );

        min_length = size_t( val );
    }

    void args_t::parse_minqscore( const char * str )
    {
        long val = atoi( str );

        if ( val < 0 )
            ERROR( "min q-score expected a non-negative integer, had: %s", str );

        min_qscore = size_t( val );
    }

    void args_t::parse_mode( const char * str )
    {
        int mode = atoi( str );

        if ( mode < 0 || mode > 7 )
            ERROR( "mode must be an integer in [0, 7], had: %s", str );

        split = ( mode & 1 );
        hpoly = ( mode & 2 );
        ambig = ( mode & 4 );
    }

    void args_t::parse_split()
    {
        split = true;
    }

    void args_t::parse_hpoly()
    {
        hpoly = true;
    }

    void args_t::parse_ambig()
    {
        ambig = true;
    }

    void args_t::parse_json()
    {
        json = true;
    }

    void args_t::parse_punch( const char * str )
    {
        const size_t len = strlen( str );

        if ( len != 1 )
            ERROR( "punch character must have a length of 1, had %s", str );

        punch = str[0];
    }

    void args_t::parse_tag( const char * str )
    {
        const int nvar = sscanf( str, "%256s", tag );

        if ( nvar != 1 )
            ERROR( "failed to process tag argument %s", str );

        tag_length = strlen( tag );
    }

    void args_t::parse_tagmismatch( const char * str )
    {
        long val = atoi( str );

        if ( val < 0 )
            ERROR( "maximum tag mismatch expected non-negative integer, had: %s", str );

        tag_mismatch = size_t( val );
    }

    void args_t::parse_format( const char * str )
    {
        if ( !strcmp( str, "FASTA" ) )
            format = FASTA;
        else if ( !strcmp( str, "FASTQ" ) )
            format = FASTQ;
        else
            ERROR( "invalid format %s", str );
    }
}
