// uses NTL
//   http://www.shoup.net/ntl

#include "BCHCode.h"

#define BCH_GF2E_DEGMAX 16

BCHCode::BCHCode(long N1, long T1)
// N1 = code length (after encodeing)
// T1 = max number of errors to be corrected
// assume GF2E has been initialized;
// Modulus of GF2E must be such that 2 is primitive
: N(N1), T(T1)
{
    long i,m(GF2E::degree()),M(1<<m);
    GF2X f;
    GF2E a(1),b;

    if(m > BCH_GF2E_DEGMAX)
        Error("GF2E degree too large");
    if(N >= M) Error("N is too large");

    E.SetLength(M);
    SetX(f);
    conv(b,f);
    for(i=0; i<N; i++) {
        E[LongFromGF2X(rep(a))] = i;
        a /= b;
    }

    set(G);
    sqr(a,b);
    for(i=0; i<T; i++) {
        MinPolyMod(f,rep(b),GF2E::modulus());
        if(!divide(G,f)) G *= f;
        b *= a;
    }
    K = N - deg(G);
    if(K<=0) Error("T is too large");

}

void BCHCode::encode(Vec<GF2>& t, const Vec<GF2>& s)
// t = encoding of s
// t.length = multiple of N; &t==&s is allowed
{
    Vec<GF2> c;
    GF2X f,g,r;
    conv(f,s);
    t.SetLength(0);
    while(!IsZero(f)) {
        trunc(g,f,K);
        reverse(g,g,K-1);
        g <<= deg(G);
        rem(r,g,G);
        g += r;
        VectorCopy(c,g,N);
        reverse(c,c);
        t.append(c);
        f >>= K;
    }
}

long BCHCode::decode(GF2X& c, const GF2X& d)
// c = decoding of d by euclidean algorithm
// return number of corrected bits if successful
// return -1 if unable to correct
// assume d.length == N; &c==&d is allowed
// c.length = K
// reference: O. Pretzel
//   "Error Correcting Codes and Finite Fields" section 16.2
{
    long i,j,M(deg(G));
    GF2X g;
    GF2E a;
    GF2EX s,t,v(1),u,q,r;
    Vec<Pair<GF2EX, long> > f;

    conv(r,d);
    s.SetLength(T<<1);
    SetX(g);
    for(i=j=0; i<T; i++, j+=2) {
        conv(a,g);
        eval(s[j],r,a);// syndrome
        sqr(s[j+1], s[i]);
        g <<= 2;
    }
    s.normalize();
    RightShift(c,d,M);
    if(IsZero(s)) return 0;
    SetCoeff(t, T<<1);
    while(deg(s) >= T) {// euclid
        DivRem(q,r,t,s);
        mul(t,q,v);
        sub(t,u,t);
        u = v;
        v = t;
        t = s;
        s = r;
    }
    if(IsZero(s) || IsOne(v)) return -1;
    MakeMonic(v);
    CanZass(f,v);// error locater
    if(f.length() < deg(v)) return -1;
    c.SetLength(K);
    for(i=0; i<f.length(); i++) {
        j = E[LongFromGF2X(rep(ConstTerm(f[i].a)))];
        if(j>=M) c[j-M]--;
    }
    return f.length();
}

long BCHCode::decode(Vec<GF2>& s, const Vec<GF2>& t)
// s = decoding of t; &s==&t is allowed
// returns 0 if successful, -1 if any error exist
{
    long k(0);
    GF2X f,g;
    Vec<GF2> d;
    conv(f,t);
    s.SetLength(0);
    while(!IsZero(f)) {
        trunc(g,f,N);
        reverse(g,g,N-1);
        if(decode(g,g) < 0) k=-1;
        VectorCopy(d,g,K);
        reverse(d,d);
        s.append(d);
        f >>= N;
    }
    return k;
}