// uses NTL
//   http://www.shoup.net/ntl

#include "RSCode.h"

#define RS_GF2E_DEGMAX 16

RSCode::RSCode(long N1, long K1)
// N1 = code length (after encodeing)
// K1 = message length (before encoding)
// assume GF2E has been initialized;
// Modulus of GF2E must be such that 2 is primitive
: N(N1), K(K1)
{
    long i,m(GF2E::degree()),L(1<<m),M(N-K);
    GF2X f;
    GF2E a(1),b;
    GF2EX g;

    if(m > RS_GF2E_DEGMAX)
        Error("GF2E degree too large");
    if(N >= L) Error("N is too large");
    if(M <= 1) Error("K is too large");

    SetX(f);
    conv(b,f);
    E.SetLength(L);
    for(i=0; i<N; i++) {
        E[LongFromGF2X(rep(a))] = i;
        a /= b;
    }

    set(G);
    SetCoeff(g,1);
    SetCoeff(g,0,b);
    for(i=0; i<M; i++) {
        G *= g;
        g[0] *= b;
    }
}

void RSCode::encode(Vec<GF2>& t, const Vec<GF2>& s)
// t = encoding of s
// t.length = multiple of N
{
    long i(K),m(GF2E::degree());
    GF2X f,g;
    GF2EX b,c,r;
    Vec<GF2> v;
    conv(f,s);
    t.SetLength(0);
    b.SetLength(K);
    while(!IsZero(f)) {
        trunc(g,f,m);
        conv(b[--i], g);
        f >>= m;
        if(!IsZero(f) && i>0) continue;
        while(i>0) clear(b[--i]);
        LeftShift(c,b,deg(G));
        rem(r,c,G);
        c += r;
        for(i=N-1; i>=0; i--) {
            conv(g,c[i]);
            VectorCopy(v,g,m);
            t.append(v);
        }
        i=K;
    }
}

long RSCode::decode(GF2EX& c, const GF2EX& d)
// c = decoding of d by euclidean algorithm
// return number of corrected chars if successful
// return -1 if unable to correct
// assume d.length == N; &c==&d is allowed
// c.length = K
// reference: O. Pretzel
//   "Error Correcting Codes and Finite Fields" section 17.7
{
    long i,j,M(deg(G)),T(M>>1);
    GF2X g;
    GF2E a,b;
    GF2EX s,t,v(1),u,q,r;
    Vec<Pair<GF2EX, long> > f;

    s.SetLength(M);
    SetX(g);
    for(i=0; i<M; i++) {
        conv(a,g);
        eval(s[i],d,a);// syndrome
        g <<= 1;
    }
    s.normalize();
    RightShift(c,d,M);
    if(IsZero(s)) return 0;
    SetCoeff(t,M);
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
    diff(t,v);
    MakeMonic(v);
    CanZass(f,v);// error locater
    if(f.length() < deg(v)) return -1;
    c.SetLength(K);
    for(i=0; i<f.length(); i++) {
        j = E[LongFromGF2X(rep(ConstTerm(f[i].a)))];
        if(j<M) continue;
        eval(a, s, ConstTerm(f[i].a));
        eval(b, t, ConstTerm(f[i].a));
        a /= b;
        c[j-M] -= a;
    }
    return f.length();
}

long RSCode::decode(Vec<GF2>& s, const Vec<GF2>& t)
// s = decoding of t; &s==&t is allowed
// returns 0 if successful, -1 if any error exist
{
    long i(N),k(0),m(GF2E::degree());
    GF2X f,g;
    GF2EX c,d;
    Vec<GF2> u;
    conv(f,t);
    s.SetLength(0);
    d.SetLength(N);
    while(!IsZero(f)) {
        trunc(g,f,m);
        conv(d[--i], g);
        f >>= m;
        if(!IsZero(f) && i>0) continue;
        while(i>0) clear(d[--i]);
        if(decode(c,d) < 0) k=-1;
        for(i=K-1; i>=0; i--) {
            conv(g,c[i]);
            VectorCopy(u,g,m);
            s.append(u);
        }
        i=N;
    }
    return k;
}