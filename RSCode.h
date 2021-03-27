// uses NTL
//   http://www.shoup.net/ntl

#ifndef __RSCode_h__
#define __RSCode_h__

#include<NTL/GF2EXFactoring.h>
using namespace NTL;

struct RSCode {// Reed-Solomon code
    long N;// code length (after encoding)
    long K;// message length (before encoding)
    GF2EX G;// generator polynomial
    Vec<long> E;// table of error location
    RSCode(long, long);
    void encode(Vec<GF2>&, const Vec<GF2>&);
    long decode(GF2EX&, const GF2EX&);
    long decode(Vec<GF2>&, const Vec<GF2>&);
};

long LongFromGF2X(const GF2X&);
void conv(Vec<GF2>&, const std::string&);
void conv(std::string&, const Vec<GF2>&);
void BuildPrimitive(GF2X&, long);
void AddNoise(Vec<GF2>&, const Vec<GF2>&, double);

#endif // __RSCode_h__