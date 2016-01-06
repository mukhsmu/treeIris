//
// ROOT analyzer
//
// CSI detector handling
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

#include "HandleV1190.h"
#include "Globals.h"

extern TFile* treeFile;

extern TTree* tree;

static uint32_t tRef=0, previousSN=0, tRF =0;
int gV1190nitems;

//const int Nchannels = 256;
//const int binlimit = 4196;
//TH1D * hV1190_T[Nchannels] = {NULL}; // T
//void HandleBOR_V1190(int run, int time);
//void HandleEOR_V1190(int run, int time);
const double timeSlope = 1; // ns/channel 
const uint32_t Nchannels   = 64+128+128+64+64+64;
const uint32_t choffset[] = {0, 64, 192, 320, 384, 448}; // The new order May 7/2015 !

void HandleV1190(TMidasEvent& event, void* ptr, ITdc* ptdc, int nitems, int MYLABEL)
//void HandleV1190(TMidasEvent& event, void* ptr, tdc_t* ptdc, int nitems, int MYLABEL)
//void HandleV1190(TMidasEvent& event, void* ptr, int nitems, int MYLABEL)
{
	ITdc times;
	uint32_t *data;
	int    i, debug = 0, debug1 = 0; 
	// int eventId = event.GetEventId();
	
	data = (uint32_t *) ptr;
	//if ( MYLABEL == 0 ) {
	//	if (hV1190_T[0] == NULL) HandleBOR_V1190(gRunNumber, 0);    
		// Reset Waveform histos before scanning ALL banks
		//    hV1190_T[channel]->Reset();
		
		// published var for HandlePHYSICS
		gV1190nitems = nitems;
		uint32_t geo, evtcnt=0, tdc=0, bunchid, evtid, channel;
		uint32_t modchannel, measure, wordcnt, evtidt, errflag, trailer;
		// Loop over all the banks
		if (nitems != 0) {
       		if (debug) printf("V1190_ Evt#:%d nitems:%d bk#%d\n", event.GetSerialNumber(), nitems, MYLABEL); 
			for (i=0 ; i<nitems ; i++) {
				if (debug) 	printf("items:%d - data[%d]: %d / 0x%x\n", Nchannels, i, data[i], data[i]);
				switch (data[i] & 0xF8000000) {
					case 0x40000000: // Global Header
						geo = data[i] & 0x1F;
						evtcnt = 0x3FFFFF & (data[i] >> 5);
						break;
					case 0x08000000:  // TDC Header
						tdc  = 0x3 & (data[i] >> 24);
						bunchid = data[i] & 0xFFF;
						evtid   = 0xFFF & (data[i] >> 12);
						break;  
					case 0x00000000:  // TDC measurement
						channel = 0x3F & (data[i] >> 19);
				 	  	channel = (channel/16)*16+15-channel%16; //Flipping channels 0->15; 1->14;...;15->0; 31->16; 47->32; 63;48.
 	  					modchannel = (tdc/2 * 64) + channel;  // changed AS May 15,2013
						// modchannel = (tdc * 32) + channel; // original
						measure = data[i] & 0x7FFFF;
						(data[i] & 0x04000000) ? trailer = 1 : trailer = 0;
						//if (modchannel<Nchannels) 	hV1190_T[modchannel]->Fill((float)measure, 1.);
 	  					// Deal with delta time in reference to the STOP (tRef)
	  					if (modchannel+choffset[MYLABEL] == 14 && event.GetSerialNumber() != previousSN) {
	    					previousSN = event.GetSerialNumber();
	    					// Extract tRef (ICT_@0 (1))
	    					tRef = measure;
	  	  				}  
	  					if ((modchannel+choffset[MYLABEL] == 12) && measure < 19500){ //temp fix to get the 3rd RF peak
	   	    				tRF = measure; //RF signal
							//	 printf("trf signal =%d",measure);
	  					}
 	    				// hTime[modchannel+choffset[MYLABEL]]->Fill( ((float)measure - (float)tRef)*timeSlope, 1.);
	   	    			// hTimeRF[modchannel+choffset[MYLABEL]]->Fill( ((float)measure - (float)tRF)*timeSlope, 1.);
	   					times.timeRaw[modchannel+choffset[MYLABEL]] = measure;
	   					times.timeRF[modchannel+choffset[MYLABEL]] = ( ((float)measure - (float)tRF)*timeSlope);
	    				// hhits->Fill( modchannel+choffset[MYLABEL], 1.);
	    				// hTall->Fill( ((float)measure - (float)tRef)*timeSlope, 1.);
	    				// hTall2D->Fill(modchannel+choffset[MYLABEL],((float)measure - (float)tRef)*timeSlope);
						if (debug1) {
							printf("evt#:%d trailer:%d bk:%d tdc:%d channel:%d modch:%d H#:%d measure:%d/0x%x tRef:%d/0x%x dT:%d\n"
							,evtcnt, trailer, MYLABEL, tdc, channel, modchannel, modchannel+choffset[MYLABEL], measure, measure, tRef, tRef, measure-tRef);
						}
						break;
					case 0x18000000: // TDC Trailer
						wordcnt = data[i] & 0xFFF;
						evtidt  = 0xFFF & (data[i] >> 12);
						break;
					case 0x20000000: // TDC Error
						errflag = data[i] & 0x7FFFF;
						break;
				} // switch
			} // for loop
		} // nitems!=0 
	*ptdc=times;	
	//tree->Fill();
}

//---------------------------------------------------------------------------------
//void HandleBOR_V1190(int run, int time)
//void HandleBOR_V1190(int run, int time, tdc_t *ptdc)
void HandleBOR_V1190(int run, int time, ITdc *ptdc)
{
	tree->Branch("tdc","ITdc",ptdc,32000,99);
	//tree->Branch("det","IDet",pdet,32000,99);
// 	char label[32], sig[32];
   	for (int i=0;i<512;i++){
     	ptdc->timeRF[i]=0.;
 	}

//  	// Booking 
//  	hV1190_T[0] = (TH1D*)gDirectory->Get("Time00");
//  	if (hV1190_T[0] == 0) {
//    	if (gOutputFile) {
// 			gOutputFile->cd();
//    		// Make a Time directory and cd to it.
//     		TDirectory* Time_dir = gOutputFile->mkdir("Time");      
//     		Time_dir->cd();
//
//      		for (int channel=0;channel<Nchannels;channel++) {
//				sprintf(label, "Time%02d", channel);
//				sprintf(sig, "Time%03d", channel);
//				hV1190_T[channel] = new TH1D(label, sig, binlimit, 0, binlimit);
//      		}
//    	}
//  	}
}

//---------------------------------------------------------------------------------
void HandleEOR_V1190(int run, int time)
{
	printf(" in V1190 EOR\n");
}
