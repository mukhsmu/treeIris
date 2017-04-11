//
// ROOT analyzer
//
// DAQ Stat handling
//
// $Id$
//

/// \mainpage
///
/// \section intro_sec Introduction
///
///
/// \section features_sec Features
///  
///   state = gOdb->odbReadInt("/runinfo/State");

#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <assert.h>
#include <signal.h>

//#include "midasServer.h"
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

// published var for HandlePHYSICS
int gSTATnitems;

const int Nchannels = 8;

void HandleBOR_STAT(int run, int time);
void HandleEOR_STAT(int run, int time);

void HandleSTAT(TMidasEvent& event, void* ptr, int nitems, int MYLABEL)
{
	uint32_t *data;
	int    i, debug = 1, debug1 = 0; 
	int eventId = event.GetEventId();
	
	data = (uint32_t *) ptr;
	if ( MYLABEL == 0 ){
		
		// publishig var for HandlePHYSICS
		gSTATnitems = nitems;
		
		// Loop over all the banks
		
		if (nitems != 0) {
			if (debug) printf(" STAT Evt#:%d nitems:%d\n", event.GetSerialNumber(), nitems); 
			if (nitems > Nchannels) return;
		} // nitems != 0
	} 
}

//typedef struct {
//  char label[16];
//  uint32_t nbin;
//  uint32_t lowbin;
//  uint32_t highbin;
//} Hdef;
//
//Hdef hdef[Nchannels] = { {"version", 1000, 0, 0xffffffff}
//		  ,{"TrigNo", 1000, 0,  0xffffffff}, {"TrigTS", 1000, 0,  0xffffffff}
//		  ,{"Rd StartT", 1000, 0,  0xffffffff}, {"Rd EndT", 1000, 0,  0xffffffff}
//		  ,{"Trig Latency", 300, 0,299}
//		  ,{"Rd eTime", 3000, 0, 2999}, {"Busy eTime", 3000, 0, 2999}
//};

//---------------------------------------------------------------------------------
void HandleBOR_STAT(int run, int time)
{
  printf(" in Stat BOR\n");
}

void HandleEOR_STAT(int run, int time)
{
  printf(" in Stat EOR\n");
}
