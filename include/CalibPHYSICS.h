// CalibPHYSICS.h

#ifndef CalibPHYSICS_H
#define CalibPHYSICS_H
#include <TObject.h>
#include <TClass.h>
#include <string>

//Extern
//extern int gPHYSICSnitems;
class CalibPHYSICS : public TObject {
	public:
		CalibPHYSICS(); 
		virtual ~CalibPHYSICS() {} //! 
		
		std::string installPath;
		std::string fileGates;
		std::string fileELoss;
		std::string filePdedx;
		std::string fileDdedx;
		std::string fileIdedx[3];
		std::string fileRunDepPar[3];
		std::string fileGeometry;
		std::string fileTCorrIC;

		Bool_t boolICGates;
		Bool_t boolGates;
		Bool_t boolELoss;
		Bool_t boolPdedx;
		Bool_t boolDdedx;
		Bool_t boolIdedx[3];
		Bool_t boolRunDepPar[3];
		Bool_t boolGeometry;
		Bool_t boolEssential;
		Bool_t boolTCorrIC;

		//virtual void ReadCalibPar(char* line);
		virtual void ReadFilenames(char* line);
		virtual void Load(TString filename);
		virtual void Print();
		virtual void Clear();
//		ClassDef(CalibPHYSICS,1)
};

#endif
// end
