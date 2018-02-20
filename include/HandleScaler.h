// HandleScaler.h

#ifndef HandleScaler_H
#define HandleScaler_H
#include "IScaler.h"

//Extern
extern int gScalernitems;

void HandleScaler(TMidasEvent& event, void* ptr, int wsize, int bank, IScaler * pscaler);
void  HandleBOR_Scaler(int run, int file, int time, IScaler * pscaler);
void  HandleEOR_Scaler(int run, int time);

#endif
// end

