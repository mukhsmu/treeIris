//#include <TLorentzRotation.h>
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "IParticle.h"
#include "TClass.h"
//#include "PhysicalConstants.h"

ClassImp(IParticle);

IParticle::IParticle()
{
  IParticle::Class()->IgnoreTObjectStreamer();
  fA = 0;
  fZ = 0;
  fT = 0;
  fE = 0;
  fPx = 0;
  fPy = 0;
  fPz = 0;
  fMass = 0;
  fTheta = 0;
  fTheta = 0;
  fPhi = 0;
  fDetected = 0;
}
//______________________________________________________________________________          
IParticle::IParticle(const IParticle &sp) : TObject(sp)
{
  // -- Copy constructor.                                                                 
  ((IParticle&)sp).Copy(*this);
}


//______________________________________________________________________________          
void IParticle::Copy(TObject &sp) const
{
  // -- Copy this method.                                                                 

  TObject::Copy((TObject&)sp);

  ((IParticle&)sp).fEulerPhi   = fEulerPhi;
  ((IParticle&)sp).fEulerTheta = fEulerTheta;
  ((IParticle&)sp).fEulerPsi   = fEulerPsi;

  //  ((IParticle&)sp).Clear(); //Why do we need this?                                
}

//_____________________________________________________________________________
Double_t IParticle::Boost(TVector3 v){//boost by TVector v
  fPVec.SetPxPyPzE(fPx,fPy,fPz,fE);
  fPVec.Boost(v);
  return Set4Vec(fPVec);
}

//______________________________________________________________________________          
Double_t IParticle::Set4Vec(TLorentzVector lv){
  //  fPVec.Setlv;
  fE = lv.E();
  fPx = lv.Px();
  fPy = lv.Py();
  fPz = lv.Pz();
  fP = sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
  //fMass = sqrt(fE*fE-fP*fP);
  fT = fE - sqrt(fE*fE-fP*fP);
  fPVec.SetPxPyPzE(fPx,fPy,fPz,fE);
  fTheta = fPVec.Theta()*TMath::RadToDeg();
  fPhi = fPVec.Phi()*TMath::RadToDeg();
  return sqrt(fE*fE-fP*fP)-fMass; //return the difference between the invariant mass and fMass (should be 0)
}

//______________________________________________________________________________          
void IParticle::Clear(Option_t *option)
{
  // --                                                                                   
  //                                                                                      

  SetA(0);
  SetZ(0);

 
  //fMass =0;
     fT  = 0;                                                                             
     fE  = 0;                                                                             
     fPx = 0;                                                                             
     fPy = 0;                                                                             
     fPz = 0;        
     fTheta = 0;
     fPhi = 0;       
     fThetacm = 0;                                                              
     fDetected = 0;
}


// //______________________________________________________________________________       
// Double_t IParticle::GetBeta()                                                      
// {                                                                                      
//   TLorentzVector vec(fE,fPx,fPy,fPz);                                                  
//   return vec.Beta();                                                                   
// }            


//______________________________________________________________________________          
TVector3 IParticle::GetBeta()
{
  // -- Get beta for this particle.                                                       

  Get4Vec();
  fBeta.SetXYZ(fPx/fE,fPy/fE,fPz/fE);
  return fBeta;
}


//______________________________________________________________________________          
TLorentzVector IParticle::Get4Vec()
{
  // -- Get the 4-vector for this particle.                                               

  fPVec.SetPxPyPzE(fPx,fPy,fPz,fE);
  return fPVec;
}


// //______________________________________________________________________________       
// Double_t IParticle::GetRelVel()                                                    
// {                                                                                      
//   //                                                                                   

//   //  T05133Event *event = (T05133Event*)fEvent.GetObject();                           
//   TVector3 cmVec = fEvent->GetVelCM();                                                 

//   Double_t sMass = fEvent->fSp.fA*kAmu;                                                
//   TVector3 sVec(fEvent->fSp.fPx/sMass,fEvent->fSp.fPy/sMass,fEvent->fSp.fPz/sMass);    

//   Double_t mass = fA*kAmu;                                                             
//   TVector3 hVec(fPx/mass,fPy/mass,fPz/mass);                                           
//   TVector3 vec = hVec-cmVec;                                                           

//   return vec.Mag();                                                                    
// }                                                                                      


//______________________________________________________________________________          
Double_t IParticle::GetTCM(TVector3 boost)
{
  // -- Get relativistic kinetic energy for this particle.                                

  Double_t tCM=0;

  // Initialize the Lorentz Transformation with the given boost.                          
  // TLorentzRotation rotCM(boost);
  // Invert the transformation as we are going to the rod frame.                          
  // rotCM.Invert();
  // Get the particle's 4-vector.                                                         
  TLorentzVector vec = Get4Vec();
  // Transform 4-vector into COM system.                                                  
  // vec = rotCM * vec;
  // In the COM system the spatial momentum is 0.                                         
  // The relativistic K.E. is then Etot - Rest Energy.                                    
  tCM = vec.E() - fMass;

  return tCM;
}

//______________________________________________________________________________          
TVector3 IParticle::GetVelVec()
{
  //                                                                                      

  //  Double_t mass = fA*kAmu;
  fVelVec.SetXYZ(fPx/fMass,fPy/fMass,fPz/fMass);

  return fVelVec;
}

//______________________________________________________________________________          
Double_t IParticle::GetVelMag()
{
  //                                                                                      

  TVector3 vec = GetVelVec();
  return vec.Mag();
}


//______________________________________________________________________________          
void IParticle::SetName(const Char_t *name)
{
  // -- Change (i.e. set) the name of the IParticle.                                  
  //                                                                                      

  fName = name;
}

//______________________________________________________________________________          
void IParticle::SetNameTitle(const Char_t *name, const Char_t *title)
{
  // -- Change (i.e. set) all the IParticle parameters (name and title).              
  //                                                                                      

  fName  = name;
  fTitle = title;
}

//______________________________________________________________________________          
void IParticle::SetTitle(const Char_t *title)
{
  // -- Change (i.e. set) the title of the IParticle.                                 
  //                                                                                      

  fTitle = title;
}


//Make it a method for a particle class.
//Double_t IParticle::eloss(Double_t ein, Double_t th , TF1 *func)//initial energy and thickness are given as arguments 
// {

//  // for (Int_t )

//  //Energy loss calculation
//  Double_t dx =0.1; //in microns
//  Double_t de = 0; //energy loss
//  Double_t en= ein; //the energy variable
//  Double_t pos = 0.;// the position variable


 
//  while (pos<= th){
//    de =  (dx*func->Eval(en))/2.; //
// en = en - de;
//    de =(dx*func->Eval(en))/2.;//energy loss in dx
//    en = en -de; // energy remaining after dx

//    pos = pos +dx;

//  }
//  return ein - en;
// }


//Make it a method for a particle class.
//Double_t IParticle::thickness(Double_t ein, Double_t efi , TF1 *func)//initial energy and final energy are given as arguments, calculates target thickness 
// {

//  // for (Int_t )

//  //Energy lossx calculation
//   Double_t dx = 0; 
//   Double_t de = 0.01; //energy loss step in MeV
//  Double_t en= ein; //the energy variable
//  Double_t th = 0.;// the thickness variable

//  //Integrate numerically 
 
//  while (en >= efi){
//    dx =  de/func->Eval(en)/2.; //
//    th = th +dx;
//    en = en-de;
//   dx =  de/func->Eval(en)/2.; //
//   th = th +dx;
//  }
//  return th;
// }


