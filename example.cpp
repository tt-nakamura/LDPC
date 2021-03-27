#include "LDPCCode.h"
#include "ConvolCode.h"
#include "BCHCode.h"
#include "RSCode.h"

main() {
    double sigma(0.6);
    long i;
    std::string s,t;
    Vec<GF2> u,v;
    Vec<double> d;
    GF2X f;

    s = "Hello world of error correcting codes.";
    conv(u,s);

    std::cout << "LDPC code ";
    LDPCCode l(512,256,3);
    l.encode(v,u);
    AddNoise(d,v,sigma);
    i = l.decode(v,d);
    conv(t,v);
    std::cout << (i ? "failed" : "succeeded") << std::endl;
    std::cout << t << std::endl;

    std::cout << "convolutional code" << std::endl;
    ConvolCode c;
    c.encode(v,u);
    AddNoise(d,v,sigma);
    c.decode(v,d);
    conv(t,v);
    std::cout << t << std::endl;
 
    BuildPrimitive(f,8);
    GF2E::init(f);

    std::cout << "BCH code ";
    BCHCode b(255,18);
    b.encode(v,u);
    AddNoise(v,v,sigma);
    i = b.decode(v,v);
    conv(t,v);
    std::cout << (i ? "failed" : "succeeded") << std::endl;
    std::cout << t << std::endl;

    std::cout << "RS code ";
    RSCode r(64,32);
    r.encode(v,u);
    AddNoise(v,v,sigma);
    i = r.decode(v,v);
    conv(t,v);
    std::cout << (i ? "failed" : "succeeded") << std::endl;
    std::cout << t << std::endl;
}