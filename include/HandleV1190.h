// HandleV1190.h

#ifndef HandleV1190_H
#define HandleV1190_H
#include "TMidasEvent.h"
#include "TEvent.h"



//Extern
extern int gV1190nitems;

class tdc_t
{
	public: 
		Double_t timeRF[512];
		int timeRaw[512];
};

//void  HandleV1190(TMidasEvent& event, void* ptr, tdc_t* timeArray, int wsize, int MYLABEL);
void  HandleV1190(TMidasEvent& event, void* ptr, ITdc* timeArray, int wsize, int MYLABEL);
//void  HandleBOR_V1190(int run, int time, tdc_t *timeArray);
void  HandleBOR_V1190(int run, int time, ITdc *timeArray);
void  HandleEOR_V1190(int run, int time);

#endif
// end
