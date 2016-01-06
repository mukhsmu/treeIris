// CalibMesytec.h

#ifndef CalibMesytec_H
#define CalibMesytec_H
#include <TObject.h>
#include <TClass.h>

//Extern
//extern int gMesytecnitems;
class CalibMesytec : public TObject {
	public:
		CalibMesytec(); 
		virtual ~CalibMesytec() {} //! 

		TString installPath;
		TString fileGeometry;
		TString fileELoss;
		TString fileIC;
		TString fileCsI1;
		TString fileCsI2;
		TString fileSd1r;
		TString fileSd1s;
		TString fileSd2r;
		TString fileSd2s;
		TString fileSur;
		TString fileSus;
		TString fileYd;
		TString fileYu;

		Bool_t boolGeometry;
		Bool_t boolELoss;
		Bool_t boolIC;
		Bool_t boolCsI1;
		Bool_t boolCsI2;
		Bool_t boolSd1r;
		Bool_t boolSd1s;
		Bool_t boolSd2r;
		Bool_t boolSd2s;
		Bool_t boolSur;
		Bool_t boolSus;
		Bool_t boolYd;
		Bool_t boolYu;
		Bool_t boolASCII;

		//virtual void ReadCalibPar(char* line);
		virtual void ReadFilenames(char* line);
		virtual void Load(TString filename);
		virtual void Print();
		virtual void Clear();
//		ClassDef(CalibMesytec,1)
};

#endif
// end
