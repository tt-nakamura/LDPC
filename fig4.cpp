#include "LDPCCode.h"
#include "ConvolCode.h"
#include "BCHCode.h"
#include<fstream>

main() {
    long N(512),K(256),W(3),IMAX(32);
    long n(12),M(100);
    double s1(0.6), s2(1.2);
    double sigma, e, e1[3], e2[3], ds((s2-s1)/n);
    long i,j,k;
    Vec<GF2> u,v;
    Vec<double> y;
    std::ofstream f("fig4.txt");
    GF2X g;

    LDPCCode l(N,K,W,IMAX);
    ConvolCode c;
    BuildPrimitive(g,4);
    GF2E::init(g);
    BCHCode b(15,2);

    for(i=0; i<=n; i++) {
        sigma = s1 + i*ds;
        std::cout << sigma << ' ';
        f << sigma << ' ';
        for(k=0; k<3; k++) e1[k] = e2[k] = 0;
        for(j=0; j<M; j++) {
            random(u,K);
            // LDPC
            l.encode(v,u);
            AddNoise(y,v,sigma);
            l.decode(v,y);
            e = weight(v-=u);
            e /= v.length();
            e1[0] += e;
            e2[0] += e*e;
            // Convol
            c.encode(v,u);
            AddNoise(y,v,sigma);
            c.decode(v,y);
            v.SetLength(u.length());
            e = weight(v-=u);
            e /= v.length();
            e1[1] += e;
            e2[1] += e*e;
            // BCH
            b.encode(v,u);
            AddNoise(v,v,sigma);
            b.decode(v,v);
            v.SetLength(u.length());
            e = weight(v-=u);
            e /= v.length();
            e1[2] += e;
            e2[2] += e*e;
        }
        for(k=0; k<3; k++) {
            e1[k] /= M;
            e2[k] = sqrt(e2[k]/M - e1[k]*e1[k]);
            std::cout << e1[k] << ' ';
            std::cout << e2[k] << ' ';
            f << e1[k] << ' ';
            f << e2[k] << ' ';
        }
        std::cout << std::endl;
        f << std::endl;
    }
}