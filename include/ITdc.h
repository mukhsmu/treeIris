// ITdc.h

#ifndef ITdc_H
#define ITdc_H
#include <TObject.h>
#include <TClass.h>

//Extern
//extern int gV1190nitems;

class ITdc : public TObject {
 	public:
  	// Int_t EventNumber;
  		ITdc(); //! Create
  		virtual ~ITdc() {} //!  
		Double_t timeRF[512];
		int timeRaw[512];

  		void Clear();
  		ClassDef(ITdc,1)
};
#endif
// end
