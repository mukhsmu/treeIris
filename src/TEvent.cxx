// Author: Alisher Sanetullaev  TRIUMF, 2012/10/01                                           
#include <TClonesArray.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <fstream>

//#include "PhysicalConstants.h"
#include "TEvent.h"
//#include "IParticle.h"
#ifndef kAmu
#define kAmu 931.494013
#endif

//ClassImp(TEvent);
//ClassImp(IParticle);
//using namespace TMath;

using namespace std;

//ofstream fout;

TEvent::TEvent() {
  //  TBranch       *b_fqV;                                            // */

   TTree        *fChain;                       // pointer to the analyzed TTree or TChain

  TEvent::Class()->IgnoreTObjectStreamer();
  //IParticle::Class()->IgnoreTObjectStreamer();   
  //fPart = new TParticle();        //

  // fEYY1 =0;             //
  // fCsI=0;              //
  // fSd1=0;              //
  // fSd2=0;              //
  // fIC=0;               //

  fLP=0; //Light particle energy //
  fHP=0; //Heavy particle energy //
  fQv=0;                         //
  fThetacm=0;                         //
  fEB=0; //Beam energy           //


}
TEvent::TEvent(Int_t k) {fRun = k;}
TEvent::TEvent(const TEvent &evt) : TObject(evt) //
{
  // -- Copy constructor.                                                                
  ((TEvent&)evt).Copy(*this);
}

Int_t TEvent::Clear(){
  // fEYY1 =0;             //
  // fCsI=0;              //
  // fSd1=0;              //
  // fSd2=0;              //
  // fIC=0;               //

  fLP=0; //Light particle energy //
  fHP=0; //Heavy particle energy //
  fQv=0;                         //
  fThetacm=0;                         //
  fEB=0; //Beam energy           //
  // fPart.Clear();
return 1;
}

Int_t TEvent::Calculate(){

  //  fEB = fHP;
  fQv = 1.;
   return 1;
}
