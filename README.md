TN93
----

This is a simple program meant to compute pairwise distances between aligned 
nucleotide sequences in sequential FASTA format using the 
[Tamura Nei 93 distance](http://www.ncbi.nlm.nih.gov/pubmed/8336541).


To build, you need to use cmake. Type 
```
git clone https://github.com/veg/tn93.git
cd tn93
cmake .
make install
```

Note, you can use :

`cmake [-DCMAKE_INSTALL_PREFIX=/install/path DEFAULT /usr/local/] ./`

to set a different install path.

If the compiler supports OpenMP, the program will be built with multithreaded
support.

USAGE
-----

	usage: tn93 [-h] [-o OUTPUT] [-t THRESHOLD] [-a AMBIGS] [-l OVERLAP][-d COUNTS_IN_NAME] [-f FORMAT] [-s SECOND_FASTA] [-b] [-c] [-q] [FASTA]

Try it from using the example file in 'data' by typing 

	tn93 -t 0.05 -o data/test.dst data/test.fas

Output (diagnostics written to stderr, histogram written to stdout so can be redirected)

Example:

    Read 8 sequences of length 1320
    Will perform 28 pairwise distance calculations
    Progress:     100% (       7 links found,          inf evals/sec)
    {
	    "Actual comparisons performed" :28,
	    "Total comparisons possible" : 28,
	    "Links found" : 7,
	    "Maximum distance" : 0.0955213,
	    "Mean distance" : 0.0644451,
	    "Histogram" : [[0.005,0],[0.01,0],[0.015,0],[0.02,0],[0.025,0],[0.03,2],[0.035,1],[0.04,0],[0.045,1],[0.05,3],[0.055,1],[0.06,2],[0.065,2],[0.07,3],[0.075,4],[0.08,3],[0.085,3],[0.09,2],[0.095,0],[0.1,1],[0.105,0],[0.11,0],[0.115,0],[0.12,0],[0.125,0],[0.13,0],[0.135,0],[0.14,0],[0.145,0],[0.15,0],[0.155,0],[0.16,0],[0.165,0],[0.17,0],[0.175,0],[0.18,0],[0.185,0],[0.19,0],[0.195,0],[0.2,0],[0.205,0],[0.21,0],[0.215,0],[0.22,0],[0.225,0],[0.23,0],[0.235,0],[0.24,0],[0.245,0],[0.25,0],[0.255,0],[0.26,0],[0.265,0],[0.27,0],[0.275,0],[0.28,0],[0.285,0],[0.29,0],[0.295,0],[0.3,0],[0.305,0],[0.31,0],[0.315,0],[0.32,0],[0.325,0],[0.33,0],[0.335,0],[0.34,0],[0.345,0],[0.35,0],[0.355,0],[0.36,0],[0.365,0],[0.37,0],[0.375,0],[0.38,0],[0.385,0],[0.39,0],[0.395,0],[0.4,0],[0.405,0],[0.41,0],[0.415,0],[0.42,0],[0.425,0],[0.43,0],[0.435,0],[0.44,0],[0.445,0],[0.45,0],[0.455,0],[0.46,0],[0.465,0],[0.47,0],[0.475,0],[0.48,0],[0.485,0],[0.49,0],[0.495,0],[0.5,0],[0.505,0],[0.51,0],[0.515,0],[0.52,0],[0.525,0],[0.53,0],[0.535,0],[0.54,0],[0.545,0],[0.55,0],[0.555,0],[0.56,0],[0.565,0],[0.57,0],[0.575,0],[0.58,0],[0.585,0],[0.59,0],[0.595,0],[0.6,0],[0.605,0],[0.61,0],[0.615,0],[0.62,0],[0.625,0],[0.63,0],[0.635,0],[0.64,0],[0.645,0],[0.65,0],[0.655,0],[0.66,0],[0.665,0],[0.67,0],[0.675,0],[0.68,0],[0.685,0],[0.69,0],[0.695,0],[0.7,0],[0.705,0],[0.71,0],[0.715,0],[0.72,0],[0.725,0],[0.73,0],[0.735,0],[0.74,0],[0.745,0],[0.75,0],[0.755,0],[0.76,0],[0.765,0],[0.77,0],[0.775,0],[0.78,0],[0.785,0],[0.79,0],[0.795,0],[0.8,0],[0.805,0],[0.81,0],[0.815,0],[0.82,0],[0.825,0],[0.83,0],[0.835,0],[0.84,0],[0.845,0],[0.85,0],[0.855,0],[0.86,0],[0.865,0],[0.87,0],[0.875,0],[0.88,0],[0.885,0],[0.89,0],[0.895,0],[0.9,0],[0.905,0],[0.91,0],[0.915,0],[0.92,0],[0.925,0],[0.93,0],[0.935,0],[0.94,0],[0.945,0],[0.95,0],[0.955,0],[0.96,0],[0.965,0],[0.97,0],[0.975,0],[0.98,0],[0.985,0],[0.99,0],[0.995,0],[1,0]]
    }

DEPENDENCIES
------------
* gcc >= 5.0.0
* cmake >= 3.0.0

ARGUMENTS
---------

    optional arguments:
      -h, --help               show this help message and exit
      -v, --version            show tn93 version 
      -o OUTPUT                direct the output to a file named OUTPUT (default=stdout)
      -t THRESHOLD             only report (count) distances below this threshold (>=0, default=0.015)
      -a AMBIGS                handle ambigous nucleotides using one of the following strategies (default=resolve)
                               resolve: resolve ambiguities to minimize distance (e.g.R matches A);
                               average: average ambiguities (e.g.R-A is 0.5 A-A and 0.5 G-A);
                               skip: do not include sites with ambiguous nucleotides in distance calculations;
                               gapmm: a gap ('-') matched to anything other than another gap is like matching an N (4-fold ambig) to it;
                               a string (e.g. RY): any ambiguity in the list is RESOLVED; any ambiguitiy NOT in the list is averaged (LIST-NOT LIST will also be averaged);
      -g FRACTION              in combination with AMBIGS, works to limit (for resolve and string options to AMBIG)
                               the maximum tolerated FRACTION of ambiguous characters; sequences whose pairwise comparisons
                               include no more than FRACTION [0,1] of sites with resolvable ambiguities will be resolved
                               while all others will be AVERAGED (default=1.0)
      -f FORMAT                controls the format of the output unless -c is set (default=csv)
                               csv: seqname1, seqname2, distance;
                               csvn: 1, 2, distance;
                               hyphy: {{d11,d12,..,d1n}...{dn1,dn2,...,dnn}}, where distances above THRESHOLD are set to 100;
      -l OVERLAP               only process pairs of sequences that overlap over at least OVERLAP nucleotides (an integer >0, default=100):
      -d COUNTS_IN_NAME        if sequence name is of the form 'somethingCOUNTS_IN_NAMEinteger' then treat the integer as a copy number
                               when computing distance histograms (a character, default=':'):
      -s SECOND_FASTA          if specified, read another FASTA file from SECOND_FASTA and perform pairwise comparison BETWEEN the files (default=NULL)
      -b                       bootstrap alignment columns before computing distances (default = false)
                               when -s is supplied, permutes the assigment of sequences to file
      -r                       if -b is specified AND -s is supplied, using -r will bootstrap across sites
                               instead of allocating sequences to 'compartments' randomly
      -c                       only count the pairs below a threshold, do not write out all the pairs 
      -m                       compute inter- and intra-population means suitable for FST caclulations
                               only applied when -s is used to provide a second file  
      -u PROBABILITY           subsample sequences with specified probability (a value between 0 and 1, default = 1.0) 
      -0                       report distances between each sequence and itself (as 0); this is useful to ensure every sequence
                               in the input file appears in the output, e.g. for network construction to contrast clustered/unclustered
      -q                       do not report progress updates and other diagnostics to stderr 
      FASTA                    read sequences to compare from this file (default=stdin)


NOTES
-----

All sequences must be aligned and have the same length.  Only IUPAC characters are recognized (e.g. no ~). Sequence names can include copy number as in 

    >seqname:10

':' can be replaced with another character using `-d`, and sequences that have no explicit copy number are assumed to be a single copy. Copy numbers 
only affect histogram and mean calculations.


	
