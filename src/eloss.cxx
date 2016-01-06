//Alisher Sanetullaev
//Energy loss calculator using Bethe formula
//-dE/dx = 4PI/(m_e c^2) (nz^2)/beta^2 (e^2/4PI(epsilon_0)^2) [ln(2m_e c^2 beta^2/I (1-beta^2))-beta^2] 
//void EneLoss(char* filename, Double_t & th, Double_t & en); //Stopping power file, target thickness, intial energy



#include <TF1.h>
#include "../include/eloss.h"



Double_t fncBethe(Double_t *x, Double_t *par) //Bethe like function for fitting
{
  Double_t val=0;
 
  //  if (par[1]!=0) val = par[0]*log(x[0]*par[1])/x[0]; //Bethe
  //  if (par[1]>0) val = par[0]*TMath::Log(x[0]/par[1])/x[0]; //Bethe
  //val = par[0]+par[1]*x[0]+ par[2]*x[0]*x[0];//quadratic function
   val = par[0]+par[1]*x[0]+ par[2]*x[0]*x[0]+ par[3]*x[0]*x[0]*x[0]+ par[4]*x[0]*x[0]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
     // val = par[0]+par[1]*x[1]+ par[2]*x[1]*x[1]+ par[3]*x[1]*x[1]*x[1]+ par[4]*x[1]*x[1
  /***Fit

  ****************************************************************************************
*/

  return val;
}

//Make it a method for a particle class.
Double_t eloss(Double_t ein, Double_t th , TGraph *func)//initial energy and thickness are given as arguments 
{
	
	// for (Int_t )
	
	//Energy loss calculation
	Double_t dx =th/100.; //in mg/cm2
	Double_t de = 0; //energy loss
	Double_t en= ein; //the energy variable
	Double_t pos = 0.;// the position variable
	
	while (pos<= th){
	  	de =  (dx*func->Eval(en))/2.; //
	   	en = en - de;
	  	de =(dx*func->Eval(en))/2.;//energy loss in dx
	  	en = en -de; // energy remaining after dx
	  	pos = pos +dx;
	}

	return ein - en;
}

//Make it a method for a particle class.
Double_t elossFi(Double_t efi, Double_t th , TGraph *func)//final energy and thickness are given as arguments 
{

 // for (Int_t )
    if (th==0) return efi;

 //Energy loss calculation
 Double_t dx =th/100.; //in 
 Double_t de = 0; //energy loss
 Double_t en= efi; //the energy variable
 Double_t pos = 0.;// the position variable
 
 while (pos<= th){
   de =  (dx*func->Eval(en))/2.; //
en = en + de;
   de =(dx*func->Eval(en))/2.;//energy loss in dx
   en = en +de; // energy remaining after dx

   pos = pos +dx;

 }
 return en-efi;
}


//Make it a method for a particle class.
Double_t thickness(Double_t ein, Double_t efi , TGraph *func)//initial energy and final energy are given as arguments, calculates target thickness 
{

 // for (Int_t )
  if (ein<=efi) return 0.0;
 //Energy loss calculation
  Double_t dx = 0; 
  Double_t de = (ein-efi)/1000.; //energy loss step in MeV
  //  Double_t de = 0.01; //energy loss step in MeV
 Double_t en= ein; //the energy variable
 Double_t th = 0.;// the thickness variable

 //Integrate numerically 
 
 while (en > efi){
   dx =  de/func->Eval(en)/2.; //
   th = th +dx;
   en = en-de;

   dx =  de/func->Eval(en)/2.; //
   th = th +dx;
   
 }
 return th;
}


// Remaining energy calculator
Double_t erem(Double_t eLoss, TGraph *func) //energy loss function for remaining energy calculation for the given energy loss and energy loss function as arguments
{
  Double_t val;
  TGraph *gr2 = new TGraph(func->GetN(), func->GetY(), func->GetX());
   // if(draw) gr->Draw("AP");
  val = gr2->Eval(eLoss)-eLoss;
if (val<0)  val = 0.;
  return val;
}
