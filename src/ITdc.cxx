//ITdc.cxx

#include "ITdc.h"
#include <math.h>
ITdc::ITdc(){
  ITdc::Clear();
}

void ITdc::Clear()
{
	
 	for(int i=0; i<32; i++){
		TICTDC[i] = 0;
		TICTime[i] = 0;
		TICTimeRF[i] = 0;
	}
		
 	for(int i=0; i<128; i++){
		TYdTDC[i] = 0;
		TYuTDC[i] = 0;
		TYdTime[i] = 0;
		TYuTime[i] = 0;
		TYdTimeRF[i] = 0;
		TYuTimeRF[i] = 0;
	}
		
 	for(int i=0; i<32; i++){
		TSurTDC[i] = 0;
		TSurTime[i] = 0;
		TSurTimeRF[i] = 0;
		TSusTDC[i] = 0;
		TSusTime[i] = 0;
		TSusTimeRF[i] = 0;
		
		TSd1rTDC[i] = 0;
		TSd1rTime[i] = 0;
		TSd1rTimeRF[i] = 0;
		TSd1sTDC[i] = 0;
		TSd1sTime[i] = 0;
		TSd1sTimeRF[i] = 0;
		
		TSd2rTDC[i] = 0;
		TSd2rTime[i] = 0;
		TSd2rTimeRF[i] = 0;
		TSd2sTDC[i] = 0;
		TSd2sTime[i] = 0;
		TSd2sTimeRF[i] = 0;
	}
}

