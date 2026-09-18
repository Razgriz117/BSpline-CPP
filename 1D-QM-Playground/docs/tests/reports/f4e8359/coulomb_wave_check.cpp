#include "tise.hpp"
#include <cstdio>
#include <cmath>
int main(){
  double R=100.0; int l=1;
  for(double E: {0.1,0.2,0.3,0.4,0.5}){ double k=std::sqrt(2*E); double eta=-1.0/k;
    auto c=tise::evaluateCoulombFunctions(l,eta,k*R);
    printf("%g %.12e %.12e %.12e %.12e %.12e\n",E,c.F,c.Fprime,c.G,c.Gprime,tise::coulombPhaseShift(l,eta));}
}
