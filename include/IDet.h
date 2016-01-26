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
		Int_t TYdMul;
		Int_t TYdChannel[128];
		Double_t TYdEnergy[128];
		Int_t TYdNo[128];
		Int_t TYdRing[128];
		Double_t TYdTheta[128];// Yd theta angle

		Int_t TCsI1Mul;
		Int_t TCsI1Channel[16];
		Double_t TCsI1Energy[16];
		Double_t TCsI1ADC[16];
		Int_t TCsI2Mul;
		Int_t TCsI2Channel[16];
		Double_t TCsI2Energy[16];
		Int_t TCsIMul;
		Int_t TCsIChannel[16];
		Double_t TCsIEnergy[16];
		Double_t TCsI2ADC[16];

		Double_t SSB;
		
		Int_t TICChannel;
		Double_t TICEnergy;
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

		Int_t TSd1rMul;
		Int_t TSd1rChannel[32];
		Double_t TSd1rEnergy[32];
		Double_t TSd1rADC[32];
		
		Int_t TSd1sMul;
		Int_t TSd1sChannel[32];
		Double_t TSd1sEnergy[32];
		Double_t TSd1sADC[32];
		
		Int_t TSd2rMul;
		Int_t TSd2rChannel[32];
		Double_t TSd2rEnergy[32];
		Double_t TSd2rADC[32];
		Double_t TSd2rEnergyCal;
		
		Int_t TSd2sMul;
		Int_t TSd2sChannel[32];
		Double_t TSd2sEnergy[32];
		Double_t TSd2sADC[32];
		
		Double_t TSdETot;
		Double_t TSdTheta[32];
		Double_t TSdThetaCM;
		Double_t TSdPhi[32];

		Double_t TQv1; // Q-value from S3
		Double_t TQv2; // Q-value from S3
		Double_t QValue1; //Q-value from YY1 and CsI1
		Double_t QValue2; //Q-value from YY1 and CsI2
		Double_t QValue; //Q-value from YY1 and CsI2
		Double_t TBE; //Beam energy
	
	void Clear();
	ClassDef(IDet,1)
};

#endif
// end
