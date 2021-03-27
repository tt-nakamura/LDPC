// uses NTL
//   http://www.shoup.net/ntl

#include "ConvolCode.h"

#define CONVOL_DEGMAX 8
#define MINUS_INF (-HUGE_VAL)

const long ConvolCode::coeff[2] = {0171, 0133};

ConvolCode::ConvolCode(const long *c, long N)
// c = array of coefficients of generator polynomials
// N = number of generator polynomials
{
    long i,j;
    g.SetLength(N);
    for(i=m=0; i<N; i++) {
        GF2XFromLong(g[i],c[i]);
        if(deg(g[i]) > m) m = deg(g[i]);
    }
    if(m > CONVOL_DEGMAX) Error("degree too large");
    C.SetDims(N,1<<(m+1));
    for(i=0; i<N; i++)
        for(j=1; j<C.NumCols(); j++)
            if(weight(j&c[i])&1) set(C[i][j]);
}

void ConvolCode::encode(Vec<GF2>& t, const Vec<GF2>& s)
// t = encoding of s
{
    long i,j,k;
    GF2X f,h;
    conv(f,s);
    t.SetLength((deg(f)+m+1)*g.length());
    clear(t);
    for(i=0; i<g.length(); i++) {
        mul(h,f,g[i]);
        for(j=0, k=i; j<=deg(h); j++, k+=g.length())
            t[k] = h[j];// interleave
    }
}

void ConvolCode::decode(Vec<GF2>& s, const Vec<double>& t)
// output: s = decoding of t by Viterbi method
// input: t = received code with gaussian noise
//          = log(prob(s==0)/prob(s==1))
// reference: S. J. Johnson
//   "Iterative Error Correction" algorithm 4.3
{
    long i,j,k,l,h;
    long L(t.length()/g.length()), M(1<<m);
    double A[M],A1[M],a,b;
    Mat<long> B;
    B.SetDims(L,M);
    A[0] = 0;
    for(j=1; j<M; j++) A[j] = MINUS_INF;
    for(i=k=0; i<L; i++, k+=g.length()) {
        for(j=0; j<M; j++) {
            h = j|M;
            a = A[j>>1];
            b = A[h>>1];
            for(l=0; l<g.length(); l++) {
                if(IsOne(C[l][j])) a -= t[k+l];
                else a += t[k+l];
                if(IsOne(C[l][h])) b -= t[k+l];
                else b += t[k+l];
            }
            if(a>b) { A1[j] = a; B[i][j] = j>>1; }
            else    { A1[j] = b; B[i][j] = h>>1; }
        }
        for(j=0; j<M; j++) A[j] = A1[j];
    }
    s.SetLength(L);
    clear(s);
    for(i=L-1, j=0; i>=0; j=B[i--][j])
        if(j&1) set(s[i]);
}