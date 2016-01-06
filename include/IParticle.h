#ifndef __IPARTICLE_H
#define __IPARTICLE_H
//#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TObject.h>
#include "TVector3.h"
#include <TROOT.h>
//#include "TF1.h"
#ifndef kAmu
#define kAmu 931.494013
#endif

using namespace std;


class IParticle : public TObject {
 private:
  TString     fName;                                   //! Object identifier
  TString     fTitle;                                  //! Object title

  UShort_t           fA;                               //!
  UShort_t           fZ;                               //!
  Double_t fMass;                                      //!
  TVector3           fVelVec;                                        //!
  TVector3           fBeta;                                          //!
  TLorentzVector     fPVec;                                          //!

  
 protected:

  // public:


 private:
  Double_t     fEulerPhi;                     //! Phi Euler angle.
  Double_t     fEulerTheta;                   //! Theata Euler angle.
  Double_t     fEulerPsi;                     //! Psi Euler angle.

 public:
  Double_t           fT;
  Double_t           fE;
  Double_t           fPx;
  Double_t           fPy;
  Double_t           fPz;
  Double_t           fP;
  Double_t fTheta; //theta (z-axis is 0 angle)
  Double_t fThetacm; //theta (z-axis is 0 angle)
  Double_t fPhi; //phi around z-axis
  Double_t fDetected;//detected or not

 public:
    IParticle(); //! Create
    //: fName(), fTitle() {}
  virtual ~IParticle() {} //!  
IParticle(const IParticle &);                          // The copy constructor.

 const Char_t  *GetName() const { return fName; } //!
 const Char_t  *GetTitle() const { return fTitle; } //!
 void           SetName(const Char_t *name); //!
 void           SetNameTitle(const Char_t *name, const Char_t *title); //!
 void           SetTitle(const Char_t *title="");   //!
  
virtual void   Clear(Option_t *option=""); //!
  void         Copy(TObject &sp) const;                          //! The copy method.

  UShort_t     GetA() {return fA;}
  UShort_t     GetZ(){return fZ;}
  Double_t     GetMass(){return fMass;}
  TVector3     GetBeta();
  Double_t     GetTCM(TVector3 boost);
  Double_t Boost(TVector3 v);
  Double_t Set4Vec(TLorentzVector lv);
  Double_t     GetVelMag();
  TVector3     GetVelVec();
  TLorentzVector Get4Vec();

  void         SetA(UShort_t A){fA=A;}
  void         SetZ(UShort_t Z){fZ=Z;}
  void         SetMass(Double_t mass){fMass=mass;}
  /*   Double_t     GetRelVel(); */

  //Double_t  eloss(Double_t, Double_t, TF1 *);

//Double_t thickness(Double_t, Double_t, TF1 *);//initial energy and final energy are given as arguments, calculates target thickness 

  void         SetEulerPhi(Double_t phi)    {fEulerPhi=phi;}
  void         SetEulerTheta(Double_t theta){fEulerTheta=theta;}
  void         SetEulerPsi(Double_t psi)    {fEulerPsi=psi;}

  ClassDef(IParticle,4)     // A Particle class.
};

#endif

