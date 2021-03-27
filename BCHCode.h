// uses NTL
//   http://www.shoup.net/ntl

#ifndef __BCHCode_h__
#define __BCHCode_h__

#include<NTL/GF2EXFactoring.h>
using namespace NTL;

struct BCHCode {// Bose-Chaudhuri-Hocquenghem code
    long N;// code length (after encoding)
    long K;// message length (before encoding)
    long T;// max number of errors to be corrected
    GF2X G;// generator polynomial
    Vec<long> E;// table of error location
    BCHCode(long, long);
    void encode(Vec<GF2>&, const Vec<GF2>&);
    long decode(GF2X&, const GF2X&);
    long decode(Vec<GF2>&, const Vec<GF2>&);
};

long LongFromGF2X(const GF2X&);
void conv(Vec<GF2>&, const std::string&);
void conv(std::string&, const Vec<GF2>&);
void BuildPrimitive(GF2X&, long);
void AddNoise(Vec<GF2>&, const Vec<GF2>&, double);

#endif // __BCHCode_h__