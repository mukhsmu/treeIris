// Eloss.h

#ifndef EneLoss_H
#define EneLoss_H

#include "TGraph.h"
//Extern

Double_t  eloss(Double_t, Double_t, TGraph *);
Double_t  erem(Double_t, TGraph *);
Double_t  elossFi(Double_t, Double_t, TGraph *);
Double_t thickness(Double_t, Double_t, TGraph *);//initial energy and final energy are given as arguments, calculates target thickness 


Double_t fncBethe(Double_t *x, Double_t *par); //Bethe like function for fitting    
#endif
// end
