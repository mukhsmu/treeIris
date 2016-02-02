// HandleMesytec.h

#ifndef HandleMesytec_H
#define HandleMesytec_H
#include "TMidasEvent.h"
#include "TEvent.h"
#include <string>

int clearDetectors();
void  HandleMesytec(TMidasEvent& event, void* ptr, int wsize, int MYLABEL, IDet* pdet, TString CalibFile);
void  HandleBOR_Mesytec(int run, int time, IDet* pdet, std::string CalibFile);
void  HandleEOR_Mesytec(int run, int time);


#endif
// end
