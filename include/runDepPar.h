// run dependent parameters

#ifndef runDepPar_H
#define runDepPar_H

#include "TString.h"
#include <string>
//#ifndef kAmu
//#define kAmu 931.494013
//#endif
//
//const double M1H= kAmu*1.+7.2889; //1H mass (Mass excess from ie.lbl.gov) -AS
//const double M2H= kAmu*2.+13.1357; //2H mass (Mass excess from ie.lbl.gov) -AS
//const double M18O= kAmu*18.+0.7815; //18O mass (Mass excess from ie.lbl.gov) -AS
//const double M40Ar= kAmu*40.-35.0398; //40Ar mass (Mass excess from ie.lbl.gov) -AS
//const double M22Ne= kAmu*22.-8.0247; //22Ne mass (Mass excess from ie.lbl.gov) -AS
//const double M20Ne= kAmu*20.-7.0419; //20Ne mass (Mass excess from ie.lbl.gov) -AS
//const double M11Li= kAmu*11.+40.7284; //11Li mass (Mass excess from ie.lbl.gov) -AS
//const double M10C= kAmu*10.+15.6986; //10C mass (Mass excess from ie.lbl.gov) -AS
//const double M10B= kAmu*10.+12.0507; //10B mass (Mass excess from ie.lbl.gov) -AS
//const double M12C= kAmu*12.; //12C mass (Mass excess from ie.lbl.gov) -AS
//const double M12Be= kAmu*12.+25.0766; //12Be mass (Mass excess from ie.lbl.gov) -AS
//const double M20Mg = kAmu*20.+ 17.5587;  // 20Mg mass (mass excess from nndc) -- JSR
//const double M20Na = kAmu*20.+ 6.8502;   // 20Na mass (mass excess from nndc) -- JSR
//
// struct beam_t {
struct runDepPar_t {
 	double EBAC; //beam energy from accelerator
 	double energy; //beam energy behind target (also initialized with value from accelerator)
  	// double mass;// mass in MeV
  	double momentum;// momentum in MeV. Initial value calculated using EBAC, has to be properly set in HandleBOR
  	double beta;// beam beta in center-of-mass system. Initial value calculated using EBAC, has to be properly set in HandleBOR
  	double gamma;// beam gamma in center-of-mass system. Initial value calculated using EBAC, has to be properly set in HandleBOR 
  	// int A; //mass number
  	// int Z; //charge
	float ICmin; // Lower limit of ionization chamber gate for incoming beam.
	float ICmax; // Upper limit of ionization chamber gate for incoming beam.
	std::string nA;
	std::string na;
	std::string nB;
	std::string nb;
  }; //runDepPar;

// struct target_t {
// 	double mass;// mass in MeV                                                                                                                     
// 	int A; //mass number                                                                                                                           
// 	int Z; //charge                                                                                                                                
// 	double thickness; //target thickness                                                                                                      
// }; //runDepPar;                                                                                                                                  
// 
// struct ejectile_t {
// 	double mass;
// 	int A;
// 	int Z;
// };

  // void setRunDepPar(TString, beam_t *, target_t *, ejectile_t *, ejectile_t *);
  void setRunDepPar(TString, runDepPar_t *);
#endif
//end
