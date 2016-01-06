// IDet.h

#ifndef IDet_H
#define IDet_H
#include <TObject.h>
#include <TClass.h>

//Extern
//extern int gMesytecnitems;
class IDet : public TObject {
 	public:
  	// Int_t EventNumber;
  	IDet(); //! Create
	 //  det_t(const det_t &); //! 
	virtual ~IDet() {} //! 
	Double_t TYdEnergy[128];
	//Double_t TYdEnergy2;
	Double_t TCsI1Energy[16];
	Double_t TCsI2Energy[16];
	Double_t TCsIEnergy[16];
	Double_t TCsI1ADC[16];
	Double_t TCsI2ADC[16];
	Double_t TICEnergy;
	Double_t TSd2rEnergy;
	Double_t TSd2sEnergy;
	Double_t TSd1rEnergy;
	Double_t TSd1sEnergy;
	Double_t TSd1rADC;
	Double_t TSd2rEnergy2;
	Double_t TSd2sEnergy2;
	Double_t TSd2rEnergyCal;
	Double_t TSd1rEnergy2;
	Double_t TSd1sEnergy2;
	Double_t SSB;
	Double_t TSdETot;
	Double_t TSdTheta;
	Double_t TSdThetaCM;
	Double_t TSdPhi;
	Double_t TYdTheta[128];// Yd theta angle
	Double_t TQv1; // Q-value from S3
	Double_t TQv2; // Q-value from S3
	Double_t QValue1; //Q-value from YY1 and CsI1
	Double_t QValue2; //Q-value from YY1 and CsI2
	Double_t QValue; //Q-value from YY1 and CsI2
	Double_t TBE; //Beam energy
	Double_t TIC00;
	Double_t TIC01;
	Double_t TIC02;
	Double_t TIC03;
	Double_t TIC04;
	Double_t TIC05;
	Double_t TIC06;
	Double_t TIC07;
	Double_t TIC08;
	Double_t TIC09;
	Double_t TIC10;
	Double_t TIC11;
	Double_t TIC12;
	Double_t TIC13;
	Double_t TIC14;
	Double_t TIC15;
//  Double_t TTime;
//  Double_t TtRef;
	Int_t TYdMul;
	Int_t TYdChannel[128];
	//Int_t    TYdChannel2;
	Int_t TYdNo[128];
	Int_t TYdRing[128];
	Int_t TCsI1Mul;
	Int_t TCsI2Mul;
	Int_t TCsI1Channel[16];
	Int_t TCsI2Channel[16];
	//Int_t TCsI1Channel2;
	//Int_t TCsI2Channel2;
	Int_t TCsIChannel[16];
	Int_t TICChannel;
	Int_t TSd2rChannel;
	Int_t TSd2sChannel;
	Int_t TSd1rChannel;
	Int_t TSd1sChannel;
	Int_t TSd2rChannel2;
	Int_t TSd2sChannel2;
	Int_t TSd1rChannel2;
	Int_t TSd1sChannel2;
	
	void Clear();
	ClassDef(IDet,1)
};

#endif
// end
