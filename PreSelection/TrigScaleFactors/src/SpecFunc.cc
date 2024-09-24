#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "TMath.h"
#include "../interface/SpecFunc.h"

namespace MyFunc {

// Square
double sqr(double x)
{
  return x*x;
}

 // D-phi
 double phi_dist(double a, double b){
   if(fabs(a - b) > 3.14159265)
   {
    return 6.2831853 - fabs(a - b);
   }
   return fabs(a - b);
 }

 // D-R
 double dr_Func(double eta1, double eta2, double phidist12){
   return sqrt((eta1 - eta2)*(eta1 - eta2) + phidist12*phidist12);
 }

 // is Tight for Jets
 bool isTight(double Jeta, double NHaF, double NEmF, double CHaF, int NumofCons, int ChargMultypl, int NeutrMultypl)
 {
  if(fabs(Jeta)<=2.4 && NHaF<0.9 && NEmF<0.9 && NumofCons>1 && CHaF>0 && ChargMultypl>0)
    return true;
  else if(fabs(Jeta)>2.4 && fabs(Jeta)<=2.7 && NHaF<0.9 && NEmF<0.9 && NumofCons>1)
    return true;
  else if(fabs(Jeta)>2.7 && fabs(Jeta)<=3.0 && NEmF<0.99 && NEmF>0.02 && NeutrMultypl>2)
    return true;
  else if(fabs(Jeta)>3 && NEmF<0.9 && NHaF>0.02 && NeutrMultypl>10)
    return true;
  else
    return false;
 }


}