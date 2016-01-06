// HandleScaler.h

#ifndef HandleScaler_H
#define HandleScaler_H

//Extern
extern int gScalernitems;

class scal_t {
 public:
  Int_t EventNumber;
  Int_t ICSum;
  Int_t ScintSum;
  Int_t AccepTrigSum;
  Int_t FreeTrigSum;

  Int_t ICRate;
  Int_t ScintRate;
  Int_t AccepTrigRate;
  Int_t FreeTrigRate;
};

void HandleScaler(TMidasEvent& event, void* ptr, int wsize, int MYLABEL);
void  HandleBOR_Scaler(int run, int time);
void  HandleEOR_Scaler(int run, int time);

#endif
// end

