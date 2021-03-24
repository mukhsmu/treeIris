// HandleV1190.h

#ifndef HandleV1190_H
#define HandleV1190_H
#include <string>
#include <vector>
#include <array>
#include "TMidasEvent.h"
#include "TEvent.h"
#include "ITdc.h"

//Extern
extern int gV1190nitems;

class tdc_t
{
	public: 
		Double_t timeRF[512];
		Double_t timeRef[512];
		int timeRaw[512];
};

void  HandleV1190(TMidasEvent& event, void* ptr, int nitems, int bank, ITdc* ptdc);
void  HandleBOR_V1190(int run, int file, std::string configFile, ITdc *timeArray);
void  HandleEOR_V1190(int run, int time);

#endif
// end
