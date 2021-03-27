// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/GF2X.h>
#include<NTL/mat_GF2.h>
using namespace NTL;

struct ConvolCode {// convolutional code
    static const long coeff[2];
    long m;// memory order
    Mat<GF2> C;// table of output bits
    Vec<GF2X> g;// generator polynomials
    ConvolCode(const long *c=coeff, long N=2);
    void encode(Vec<GF2>&, const Vec<GF2>&);
    void decode(Vec<GF2>&, const Vec<double>&);
};

void GF2XFromLong(GF2X&, long);
long LongFromGF2X(const GF2X&);
void conv(Vec<GF2>&, const std::string&);
void conv(std::string&, const Vec<GF2>&);
void AddNoise(Vec<double>&, const Vec<GF2>&, double);