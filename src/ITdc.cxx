//ITdc.cxx
//#include "TObject.h"
//#include "TClass.h"
#include "ITdc.h"
#include <math.h>
//ClassImp(ITdc);
ITdc::ITdc(){
  ITdc::Clear();
}

void ITdc::Clear()
{
 	for(int i=0; i<512; i++){
		timeRF[i]=0.;
		timeRaw[i]=0;
	}
}

