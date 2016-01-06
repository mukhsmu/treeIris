// HandlePHYSICS.h
#ifndef HandlePHYSICS_H
#define HandlePHYSICS_H
//void  qvalue(IDet  *pdet);
void HandlePHYSICS(IDet *pdet);
void  HandleBOR_PHYSICS(int run, int time, IDet *pdet, TString CalibFile);
void  HandleEOR_PHYSICS(int run, int time);
Double_t Qvalue (Double_t ma, Double_t mA, Double_t mb,Double_t E1, Double_t Eb, Double_t thetaR);
#endif
// end

