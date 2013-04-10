#include "stringBuffer.h"
#include <stdlib.h>
#include <string.h>
#include <cfloat>
#include <iostream>

#define SWAP(A,B,C) ({(C)=(A); (A)=(B); (B)=(C);})

using namespace std;

long StringBuffer::sbDefaultLength = 16,
                   StringBuffer::sbDefaultBoost  = 16,
                                 Vector::vDefaultLength		   = 16,
                                             Vector::vDefaultBoost		   = 16,
                                                          VectorDouble::vDefaultLength		   = 16,
                                                                              VectorDouble::vDefaultBoost		   = 16;

/*---------------------------------------------------------------------------------------------------- */

StringBuffer::StringBuffer (void)
{
    sLength  = 0;
    saLength = StringBuffer::sbDefaultLength;
    sData = (char*)malloc (sizeof(char)*(saLength+1));
    sData[0] = 0;
}

/*---------------------------------------------------------------------------------------------------- */

void StringBuffer::swap (StringBuffer& src)
{
    long t;
    SWAP (sLength, src.sLength, t);
    SWAP (saLength, src.saLength, t);
    char * tc;
    SWAP (sData, src.sData, tc);
}

/*---------------------------------------------------------------------------------------------------- */

StringBuffer::~StringBuffer (void)
{
    free (sData);
}

/*---------------------------------------------------------------------------------------------------- */

void StringBuffer::appendChar (const char c)
{
    long addThis;
    if (saLength == sLength)
    {
        addThis = saLength / 8;
        if (StringBuffer::sbDefaultBoost > addThis)
            addThis = StringBuffer::sbDefaultBoost;
        saLength += addThis;
        sData = (char*)realloc (sData,sizeof(char)*(saLength+1));
    }
    sData[sLength] = c;
    sData[++sLength] = 0;
}

/*---------------------------------------------------------------------------------------------------- */

void StringBuffer::appendBuffer (const char *buffer, const long length)
{
    long addThis,
         pl = length > 0? length : strlen(buffer);

    if (pl>0)
    {
        if (saLength < sLength + pl)
        {
            addThis = saLength / 8;

            if (StringBuffer::sbDefaultBoost > addThis)
                addThis = StringBuffer::sbDefaultBoost;
            if (addThis < pl)
                addThis = pl;

            saLength += addThis;
            sData = (char*)realloc (sData,sizeof(char)*(saLength+1));
        }
        for (addThis = 0; addThis < pl; addThis++)
            sData[sLength++] = buffer[addThis];

        sData[sLength] = 0;
    }
}

/*---------------------------------------------------------------------------------------------------- */

void StringBuffer::resetString (void)
{
    sLength        = 0;
    sData[sLength] = 0;
}

/*---------------------------------------------------------------------------------------------------- */

Vector::Vector (void)
{
    vLength  = 0;
    vaLength = Vector::vDefaultLength;
    vData = (long*)malloc (sizeof(long)*vaLength);
}


/*---------------------------------------------------------------------------------------------------- */

void Vector::swap (Vector& src)
{
    long t;
    SWAP (vLength, src.vLength, t);
    SWAP (vaLength, src.vaLength, t);
    long * tc;
    SWAP (vData, src.vData, tc);
}

/*---------------------------------------------------------------------------------------------------- */

Vector::~Vector (void)
{
    free (vData);
}

/*---------------------------------------------------------------------------------------------------- */

void Vector::appendValue (const long l)
{
    long addThis;
    if (vLength == vaLength)
    {
        addThis = vaLength / 8;
        if (Vector::vDefaultBoost > addThis)
            addThis = Vector::vDefaultBoost;
        vaLength += addThis;
        vData = (long*)realloc (vData,sizeof(long)*vaLength);
    }
    vData[vLength++] = l;
}

/*---------------------------------------------------------------------------------------------------- */

void Vector::storeValue (const long v, const unsigned long l)
{
    long addThis;
    if (l >= vaLength)
    {
        addThis = l-vaLength+1;
        if (Vector::vDefaultBoost > addThis)
            addThis = Vector::vDefaultBoost;
        vaLength += addThis;
        vData = (long*)realloc (vData,sizeof(long)*vaLength);
        vLength = l+1;
    }
    vData[l] = v;
}


/*---------------------------------------------------------------------------------------------------- */

void Vector::storeVector (const Vector& v, const unsigned long l)
{
    if (l < vLength && vData[l])
        delete (Vector*)(vData[l]);

    storeValue ((long)&v, l);
}

/*---------------------------------------------------------------------------------------------------- */

void Vector::remove (const unsigned long l)
{
    if (l < vLength) {
        for (unsigned long k = l+1; k < vLength; k++) {
            vData[k-1] = vData[k];
        }
        vLength--;
    }
}


/*---------------------------------------------------------------------------------------------------- */

long Vector::extractMin (VectorDouble& values) {
    double current_min = DBL_MAX;
    long   best_index  = -1;

    for (unsigned long i = 0; i < vLength; i++) {
        double try_value = values.value(vData[i]);
        if (try_value < current_min) {
            best_index  = i;
            current_min = try_value;
        }
    }

    if (best_index >= 0) {
        current_min = vData[best_index];
        remove (best_index);
        return current_min;
    }

    return best_index;
}

/*---------------------------------------------------------------------------------------------------- */

int long_comp(const void * v1, const void *v2)
{
    long l1 = *(long*)v1,
         l2 = *(long*)v2;
    if (l1<l2)
        return -1;
    if (l1>l2)
        return 1;
    return 0;
}

/*---------------------------------------------------------------------------------------------------- */

void Vector::sort (void)
{
    qsort((void*)vData, vLength, sizeof (long), long_comp);
}

/*---------------------------------------------------------------------------------------------------- */

void Vector::resetVector (void)
{
    vLength        = 0;
}

/*---------------------------------------------------------------------------------------------------- */

VectorDouble::VectorDouble (void)
{
    vLength  = 0;
    vaLength = Vector::vDefaultLength;
    vData = (double*)malloc (sizeof(double)*vaLength);
}

/*---------------------------------------------------------------------------------------------------- */

VectorDouble::~VectorDouble (void)
{
    free (vData);
}

/*---------------------------------------------------------------------------------------------------- */

void VectorDouble::appendValue (const double l)
{
    long addThis;
    if (vLength == vaLength)
    {
        addThis = vaLength / 8;
        if (VectorDouble::vDefaultBoost > addThis)
            addThis = VectorDouble::vDefaultBoost;
        vaLength += addThis;
        vData = (double*)realloc (vData,sizeof(double)*vaLength);
    }
    vData[vLength++] = l;
}

/*---------------------------------------------------------------------------------------------------- */

void VectorDouble::storeValue (const double v, const unsigned long l)
{
    long addThis;
    if (l >= vaLength)
    {
        addThis = l-vaLength+1;
        if (VectorDouble::vDefaultBoost > addThis)
            addThis = VectorDouble::vDefaultBoost;
        vaLength += addThis;
        vData = (double*)realloc (vData,sizeof(double)*vaLength);
        vLength = l+1;
    }
    vData[l] = v;
}





