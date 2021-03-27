#include "LDPCCode.h"
#include<fstream>

main() {
    long W[] = {3,5,7};
    long C(3),N(512),K(256),n(12),M(100),IMAX(32);
    double s1(0.6), s2(1.2);
    double sigma, e, e1, e2, ds((s2-s1)/n);
    long i,j,k;
    Vec<GF2> u,v;
    Vec<double> y;
    std::ofstream f("fig3.txt");

    LDPCCode *code[C];
    for(k=0; k<C; k++)
        code[k] = new LDPCCode(N,K,W[k],IMAX);
    for(i=0; i<=n; i++) {
        sigma = s1 + i*ds;
        std::cout << sigma << ' ';
        f << sigma << ' ';
        for(k=0; k<C; k++) {
            e1 = e2 = 0;
            for(j=0; j<M; j++) {
                random(u,K);
                code[k]->encode(v,u);
                AddNoise(y,v,sigma);
                code[k]->decode(v,y);
                e = weight(v-=u);
                e /= v.length();
                e1 += e;
                e2 += e*e;
            }
            e1 /= M;
            e2 = sqrt(e2/M - e1*e1);
            std::cout << e1 << ' ';
            std::cout << e2 << ' ';
            f << e1 << ' ';
            f << e2 << ' ';
        }
        std::cout << std::endl;
        f << std::endl;
    }
}