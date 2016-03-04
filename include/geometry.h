// experiment geometry
#ifndef geometry_H
#define geometry_H
#include <TObject.h>
#include <TClass.h>
#include <string>

class geometry : public TObject {
	public:
		geometry(); 
		virtual ~geometry() {} //! 
		
		Double_t YdDistance;
		Double_t Yd1r;
		Double_t Yd2r;
		Double_t Sd1Distance;
		Double_t Sd2Distance;
		Double_t Sdr1;
		Double_t Sdr2;
		Double_t xShift;
		Double_t yShift;
		Double_t TThickness;

		void ReadGeometry(std::string);
		void Print();
};
#endif
// end
