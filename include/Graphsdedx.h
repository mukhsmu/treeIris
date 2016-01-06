// Graphsdedx.h

#ifndef GraphsH
#define GraphsH
#include <TObject.h>
#include <TClass.h>

//Extern
//extern int gMesytecnitems;
class Graphsdedx : public TObject {
	public:
		Graphsdedx(); 
		virtual ~Graphsdedx() {} //! 

		TString Ag;
		TString Al;
		TString B;
		TString D2;
		TString H2;
		TString My;
		TString P;
		TString Si;
		TString SiO2;

		Bool_t boolAg;
		Bool_t boolAl;
		Bool_t boolB;
		Bool_t boolD2;
		Bool_t boolH2;
		Bool_t boolMy;
		Bool_t boolP;
		Bool_t boolSi;
		Bool_t boolSiO2;

		//virtual void ReadCalibPar(char* line);
		virtual void ReadGraphnames(char* line);
		virtual void Load(TString filename);
		virtual void Print();
		virtual void Clear();
//		ClassDef(Graphsdedx,1)
};

#endif
// end
