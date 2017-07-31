// HandleScaler.h

#ifndef HandleScaler_H
#define HandleScaler_H
#include "IScaler.h"

//Extern
extern int gScalernitems;

void HandleScaler(TMidasEvent& event, void* ptr, IScaler * pscaler, int wsize, int MYLABEL);
void  HandleBOR_Scaler(int run, int file, int time, IScaler * pscaler);
void  HandleEOR_Scaler(int run, int time);

#endif
// end

