// ITdc.h

#ifndef ITdc_H
#define ITdc_H
#include <TObject.h>
#include <TClass.h>

class ITdc : public TObject {
 	public:
  		ITdc(); //! Create
  		virtual ~ITdc() {} //!  
		
		Int_t TICTDC[32];
		Double_t TICTime[32];
		Double_t TICTimeRF[32];
		
		Int_t TYdTDC[128];
		Double_t TYdTime[128];
		Double_t TYdTimeRF[128];
		
		Int_t TYuTDC[128];
		Double_t TYuTime[128];
		Double_t TYuTimeRF[128];
		
		Int_t TSurTDC[32];
		Double_t TSurTime[32];
		Double_t TSurTimeRF[32];
		
		Int_t TSusTDC[32];
		Double_t TSusTime[32];
		Double_t TSusTimeRF[32];
		
		Int_t TSd1rTDC[32];
		Double_t TSd1rTime[32];
		Double_t TSd1rTimeRF[32];
		
		Int_t TSd1sTDC[32];
		Double_t TSd1sTime[32];
		Double_t TSd1sTimeRF[32];
		
		Int_t TSd2rTDC[32];
		Double_t TSd2rTime[32];
		Double_t TSd2rTimeRF[32];
		
		Int_t TSd2sTDC[32];
		Double_t TSd2sTime[32];
		Double_t TSd2sTimeRF[32];

  		void Clear();
  		ClassDef(ITdc,1)
};
#endif
// end
