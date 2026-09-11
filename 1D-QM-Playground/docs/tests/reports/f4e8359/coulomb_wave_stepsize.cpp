#include "tise.hpp"
#include <cstdio>
#include <cmath>
#include <chrono>
int main(){
  double R=100.0; int l=1; double E=0.1; double k=std::sqrt(2*E); double eta=-1.0/k;
  for(double far: {50.0,100.0}) for(double h: {0.1,0.05,0.025,0.0125}){
    auto t0=std::chrono::steady_clock::now();
    auto c=tise::evaluateCoulombFunctions(l,eta,k*R,far,h);
    double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t0).count();
    printf("far=%g h=%g F=%.10e G=%.10e W=%.8f  %.2f ms\n",far,h,c.F,c.G,c.F*c.Gprime-c.Fprime*c.G,ms);}
}
