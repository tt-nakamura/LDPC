#include "LDPCCode.h"
#include<fstream>

main() {
    long N(512),K(256);
    long W(3),n(12),C(20),M(100),IMAX(32);
    double s1(0.4), s2(1), ds((s2-s1)/n);
    double sigma, c, c1, c2;
    long i,j,k;
    Vec<GF2> u,v;
    Vec<double> y;
    std::ofstream f("fig1.txt");

    LDPCCode *code[C];
    for(k=0; k<C; k++)
        code[k] = new LDPCCode(N,K,W,IMAX);
    for(i=0; i<=n; i++) {
        sigma = s1 + i*ds;
        c1 = c2 = 0;
        for(k=0; k<C; k++) {
            for(c=j=0; j<M; j++) {
                random(u,K);
                code[k]->encode(v,u);
                AddNoise(y,v,sigma);
                c += code[k]->decode_(v, y.data());
            }
            c /= M;
            c1 += c;
            c2 += c*c;
        }
        c1 /= C;
        c2 = sqrt(c2/C - c1*c1);
        std::cout << sigma << ' ';
        std::cout << c1 << ' ';
        std::cout << c2 << std::endl;
        f << sigma << ' ';
        f << c1 << ' ';
        f << c2 << std::endl;
    }
}