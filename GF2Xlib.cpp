// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/GF2X.h>
using namespace NTL;

void GF2XFromLong(GF2X& x, long a) {
    clear(x);
    x.SetLength(NumBits(a));
    for(long i=0; a; a>>=1, i++) if(a&1) set(x[i]);
}

long LongFromGF2X(const GF2X& a) {
    if(IsZero(a)) return 0;
    long i,x(1);
    for(i=deg(a)-1; i>=0; i--) {
        x<<=1;
        if(IsOne(a[i])) x|=1;
    }
    return x;
}

void conv(GF2X& x, const std::string& a) {
    GF2XFromBytes(x, (const unsigned char*)a.c_str(), a.length());
}

void conv(std::string& x, const GF2X& a) {
    x.resize(NumBytes(a));
    BytesFromGF2X((unsigned char*)x.data(), a, x.length());
}

void conv(Vec<GF2>& x, const std::string& a) {
    GF2X f;
    GF2XFromBytes(f, (const unsigned char*)a.c_str(), a.length());
    conv(x,f);
}

void conv(std::string& x, const Vec<GF2>& a) {
    GF2X f;
    conv(f,a);
    x.resize(NumBytes(f));
    BytesFromGF2X((unsigned char*)x.data(), f, x.length());
}

void BuildPrimitive(GF2X& f, long m)
// f = primitive polynomical of degree m
//     used as GF2E modulus for BCH code
// reference: Justesen and Hoholdt
//   "A Course in Error Correcting Codes" Appendix C
{
    long k(1<<m);
    if(m<2 || m>16) Error("bad degree in BuildPrimitive");
    else if(m<=4 || m==6 || m==15) k|=03;
    else if(m==5 || m==11) k|=05;
    else if(m==7 || m==10) k|=011;
    else if(m==8)  k|=035;
    else if(m==9)  k|=021;
    else if(m==12) k|=0123;
    else if(m==13) k|=033;
    else if(m==14) k|=02103;
    else if(m==16) k|=010013;
    GF2XFromLong(f,k);
}