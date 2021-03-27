#include "LDPCCode.h"
#include "ConvolCode.h"
#include "BCHCode.h"

main() {
    long K[] = {128, 256, 512, 1024};
    long N[] = {256, 512, 1024, 2048};
    long W(3),IMAX(32),M(100),C(4);
    double sigma(0.5),t,t1[3],t2[3];
    long i,j,k;
    Vec<GF2> u,v;
    Vec<double> y;
    GF2X g;

    ConvolCode c;
    BuildPrimitive(g,4);
    GF2E::init(g);
    BCHCode b(15,2);
    for(k=0; k<C; k++) {
        LDPCCode l(N[k],K[k],W,IMAX);
        for(i=0; i<3; i++) t1[i] = t2[i] = 0;
        for(j=0; j<M; j++) {
            random(u,K[k]);
            // LDPC
            t = GetTime();
            l.encode(v,u);
            t1[0] += GetTime() - t;
            AddNoise(y,v,sigma);
            t = GetTime();
            l.decode(v,y);
            t2[0] += GetTime() - t;
            // Convol
            t = GetTime();
            c.encode(v,u);
            t1[1] += GetTime() - t;
            AddNoise(y,v,sigma);
            t = GetTime();
            c.decode(v,y);
            t2[1] += GetTime() - t;
            // BCH
            t = GetTime();
            b.encode(v,u);
            t1[2] += GetTime() - t;
            AddNoise(v,v,sigma);
            t = GetTime();
            b.decode(v,v);
            t2[2] += GetTime() - t;
        }
        for(i=0; i<3; i++) {
            std::cout << t1[i]/M << ' ';
            std::cout << t2[i]/M << ' ';
        }
        std::cout << std::endl;
    }
}