#ifndef 	__STRINGBUFFER__
#define 	__STRINGBUFFER__

//__________________________________________________________________________________________

class	StringBuffer {

    char*	sData;
    unsigned long
    sLength,
    saLength;

public:

    StringBuffer 			(void);
    ~StringBuffer			(void);

    char*					getString		(void) const {
        return sData;
    }
    void					appendChar		(const char);
    void					appendBuffer	(const char*, const long = -1);
    void					resetString		(void);
    void                    swap            (StringBuffer&);
    unsigned long			length			(void) {
        return sLength;
    }
    char          getChar      (const long i) const {return sData[i];}


    static	long				sbDefaultLength,
                                sbDefaultBoost;

};


//__________________________________________________________________________________________

class	VectorDouble {

    double*	vData;

    unsigned long
    vLength,
    vaLength;

public:

    VectorDouble 			(void);
    ~VectorDouble			(void);

    void					appendValue		(const double);
    void					storeValue		(const double, const unsigned long);
    double					value			(const long idx) {
        return vData[idx];
    }
    unsigned long			length			(void) {
        return vLength;
    }


    static	long	vDefaultLength,
               vDefaultBoost;

};

//__________________________________________________________________________________________


class	Vector {

    long*	vData;

    unsigned long
    vLength,
    vaLength;

public:

    Vector 			(void);
    ~Vector			(void);

    void					appendValue		(const long);
    long                    extractMin      (VectorDouble&);
    void					resetVector		(void);
    void					remove          (const unsigned long);
    void					storeValue		(const long, const unsigned long);
    void					storeVector		(const Vector&, const unsigned long);
    void					sort			(void);
    void                    swap            (Vector&);
    long					value			(const long idx) const {
        return vData[idx];
    }
    unsigned long			length			(void) const {
        return vLength;
    }


    static	long				vDefaultLength,
                        vDefaultBoost;

};


#endif