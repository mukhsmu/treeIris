 // Author: Alisher Sanetullaev  TRIUMF, 2012/10/01                                          

#ifndef __TEVENT_H
#define __TEVENT_H

#include <TTree.h>
#include <TEventList.h>
//#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TObject.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TVector3.h>
#include "IParticle.h"
#include "IDet.h"
#include "ITdc.h"

using namespace std;

class TEvent : public TObject {
	public:
		TEvent(); //! Create
		TEvent(Int_t); //!
		TEvent(const TEvent &); //!
		virtual ~TEvent() {} //!
		
		TTree *fChain;          //! pointer to the analyzed TTree or TChain
		
		// IParticle fPart;        //
		
		// Double_t fEYY1;             //
		// Double_t fCsI;              //
		// Double_t fSd1;              //
		// Double_t fSd2;              //
		// Double_t fIC;               //
		Double_t fThetacm;          //
		
		Double_t fLP; //Light particle energy //
		Double_t fHP; //Heavy particle energy //
		Double_t fQv; //Q-value                   //
		Double_t fEB; //Beam energy           //
		Int_t Calculate(); //! Calculate derived variables
		Int_t Clear(); //! Clear values
		
	private:
		Int_t fRun; // Current run number
		ClassDef(TEvent,4) // Version
};
#endif
