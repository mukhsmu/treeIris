//
// ROOT analyzer
//
// DAQ Scaler handling
//
// $Id$
//
#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <assert.h>
#include <signal.h>

#include "TMidasEvent.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TFolder.h>
#include <TFile.h>
#include <TMath.h>

#include "Globals.h"

#include "HandleScaler.h"
#include "IScaler.h"
// published var for HandlePHYSICS
int gScalernitems;
Int_t runNo;

extern TTree *tree;

//scal_t scaler;
const int Nchannels = 32;

void HandleScaler(TMidasEvent& event, void* ptr, IScaler *pscaler, int nitems, int MYLABEL)
{
	IScaler scaler;
	uint32_t *data;
	int    i, debug = 0, debug1 = 0; 
	int eventId = event.GetEventId();
	scaler.Clear();	
	data = (uint32_t *) ptr;

    
    // publishing var for HandlePHYSICS
    gScalernitems = nitems;

    // Loop over all the banks

	if (nitems != 0) {
      	if (debug) printf(" Scaler EvtId: %d Evt#:%d nitems:%d\n", eventId, event.GetSerialNumber(), nitems); 
		if (nitems > Nchannels) return;
	  	if (MYLABEL==0){	  
      		for (i=0 ; i<nitems ; i++) {
				if (debug1) printf("items:%d - data[%d]: %d / 0x%x\n", Nchannels, i, data[i], data[i]);
	
      		} //for
	  	} //MYLABEL==0

	  	if (MYLABEL==1){	  
      		for (i=0 ; i<nitems ; i++) {
				if (debug1) printf("items:%d - data[%d]: %d / 0x%x\n", Nchannels, i, data[i], data[i]);
				scaler.ICRate = data[30];
	  			scaler.ScintRate = data[29];
	  			scaler.AccepTrigRate = data[26];
	  			scaler.FreeTrigRate = data[31];
      		} //for
	  	} //MYLABEL==0

	  	if (MYLABEL==2){	  
      		for (i=0 ; i<nitems ; i++) {
				if (debug1) printf("items:%d - data[%d]: %d / 0x%x\n", Nchannels, i, data[i], data[i]);
	  			scaler.ICSum = data[30];
	  			scaler.ScintSum = data[29];
	  			scaler.AccepTrigSum = data[26];
	  			scaler.FreeTrigSum = data[31];
      		} //for
	  	} //MYLABEL==0
	} // nitems != 0
     
    scaler.EventID = eventId;
 	
	*pscaler = scaler;
}

//---------------------------------------------------------------------------------
void HandleBOR_Scaler(int run, int gFileNumber, int time, IScaler *pscaler)
{
  	printf(" in Scaler BOR\n");
	if(gFileNumber==0) tree->Branch("scaler","Iscaler",pscaler,32000,99);
	else tree->SetBranchAddress("scaler",&pscaler);
	printf("Finished HandleBOR_Scaler\n");
}


void HandleEOR_Scaler(int run, int time)
{
  	printf(" in Scaler EOR\n");
}
