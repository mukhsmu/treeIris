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
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include "HandleV1190.h"
#include "Globals.h"

extern TFile* treeFile;

extern TTree* tree;

static uint32_t tRef=0, previousSN=0, tRF =0;
int gV1190nitems;

const double timeSlope = 1; // ns/vchannel 
const uint32_t Nvchannels   = 64+128+128+64+64+64;
const uint32_t choffset[] = {0, 64, 192, 320, 384, 448}; // The new order May 7/2015 !

Double_t timeRaw[512];
Double_t timeRef[512];
Double_t timeRF[512];

uint32_t geo, evtcnt=0, tdc=0, bunchid, evtid, vchannel;
uint32_t modchannel, measure, wordcnt, evtidt, errflag, trailer;

void HandleV1190(TMidasEvent& event, void* ptr, int nitems, int bank, ITdc* ptdc)
{
	ITdc times;
	uint32_t *data;
	int    i, debug = 0, debug1 = 0; 
	// int eventId = event.GetEventId();
	
	data = (uint32_t *) ptr;
		
	// published var for HandlePHYSICS
	gV1190nitems = nitems;
	// Loop over all the banks
	if (nitems != 0) {
   		if (debug) printf("V1190_ Evt#:%d nitems:%d bk#%d\n", event.GetSerialNumber(), nitems, bank); 
		for (i=0 ; i<nitems ; i++) {
			if (debug) 	printf("items:%d - data[%d]: %d / 0x%x\n", Nvchannels, i, data[i], data[i]);
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
					vchannel = 0x3F & (data[i] >> 19);
			 	  	vchannel = (vchannel/16)*16+15-vchannel%16; //Flipping vchannels 0->15; 1->14;...;15->0; 31->16; 47->32; 63;48.
  					modchannel = (tdc/2 * 64) + vchannel;  // changed AS May 15,2013
					// modchannel = (tdc * 32) + vchannel; // original
					measure = data[i] & 0x7FFFF;
					(data[i] & 0x04000000) ? trailer = 1 : trailer = 0;
  					// Deal with delta time in reference to the STOP (tRef)
  					if (modchannel+choffset[bank] == 14 && event.GetSerialNumber() != previousSN) {
    					previousSN = event.GetSerialNumber();
    					// Extract tRef (ICT_@0 (1))
    					tRef = measure;
  	  				}  
  					if ((modchannel+choffset[bank] == 12) && measure < 19500){ //temp fix to get the 3rd RF peak
   	    				tRF = measure; //RF signal
  					}
   					timeRaw[modchannel+choffset[bank]] = (Int_t)measure;
   					timeRef[modchannel+choffset[bank]] = ( ((Double_t)measure - (Double_t)tRef)*timeSlope);
   					timeRF[modchannel+choffset[bank]] = ( ((Double_t)measure - (Double_t)tRF)*timeSlope);
					if (debug1) {
						printf("evt#:%d trailer:%d bk:%d tdc:%d vchannel:%d modch:%d H#:%d measure:%d/0x%x tRef:%d/0x%x dT:%d\n"
						,evtcnt, trailer, bank, tdc, vchannel, modchannel, modchannel+choffset[bank], measure, measure, tRef, tRef, measure-tRef);
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

	times.Clear();
	if(bank==5){ 	// check for last bank 
		for(i=32; i<64; i++){
			if(timeRef[i]>0.){
				//times.TICTDC.push_back(timeRaw[i]);
				times.TICTime.push_back(timeRef[i]);
				times.TICTimeRF.push_back(timeRF[i]);
				times.TICTChannel.push_back(i-32);
				times.TICTMul++;
			}
		}
				
		for(i=64; i<192; i++){
			if(timeRef[i]>0.){
			//	times.TYdTDC.push_back(timeRaw[i]);
				times.TYdTime.push_back(timeRef[i]);
				times.TYdTimeRF.push_back(timeRF[i]);
				times.TYdTChannel.push_back(i-64);
				times.TYdTMul++;
			}
		}

		for(i=192; i<320; i++){
			if(timeRef[i]>0.){
			//	times.TYuTDC.push_back(timeRaw[i]);
				times.TYuTime.push_back(timeRef[i]);
				times.TYuTimeRF.push_back(timeRF[i]);
				times.TYuTChannel.push_back(i-192);
				times.TYuTMul++;
			}
		}
			
		for(i=320; i<352; i++){
			if(timeRef[i]>0.){
			//	times.TSurTDC.push_back(timeRaw[i]);
				times.TSurTime.push_back(timeRef[i]);
				times.TSurTimeRF.push_back(timeRF[i]);
				times.TSurTChannel.push_back(i-320);
				times.TSurTMul++;
			}
		}
		for(i=352; i<384; i++){
			if(timeRef[i]>0.){
			//	times.TSusTDC.push_back(timeRaw[i]);
				times.TSusTime.push_back(timeRef[i]);
				times.TSusTimeRF.push_back(timeRF[i]);
				times.TSusTChannel.push_back(i-352);
				times.TSusTMul++;
			}
		}
			
		for(i=384; i<416; i++){
			if(timeRef[i]>0.){
			//	times.TSd1rTDC.push_back(timeRaw[i]);
				times.TSd1rTime.push_back(timeRef[i]);
				times.TSd1rTimeRF.push_back(timeRF[i]);
				times.TSd1rTChannel.push_back(i-384);
				times.TSd1rTMul++;
			}
		}
		for(i=416; i<448; i++){
			if(timeRef[i]>0.){
			//	times.TSd1sTDC.push_back(timeRaw[i]);
				times.TSd1sTime.push_back(timeRef[i]);
				times.TSd1sTimeRF.push_back(timeRF[i]);
				times.TSd1sTChannel.push_back(i-416);
				times.TSd1sTMul++;
			}
		}
			
		for(i=448; i<480; i++){
			if(timeRef[i]>0.){
			//	times.TSd2rTDC.push_back(timeRaw[i]);
				times.TSd2rTime.push_back(timeRef[i]);
				times.TSd2rTimeRF.push_back(timeRF[i]);
				times.TSd2rTChannel.push_back(i-448);
				times.TSd2rTMul++;
			}
		}
		for(i=480; i<512; i++){
			if(timeRef[i]>0.){
			//	times.TSd2sTDC.push_back(timeRaw[i]);
				times.TSd2sTime.push_back(timeRef[i]);
				times.TSd2sTimeRF.push_back(timeRF[i]);
				times.TSd2sTChannel.push_back(i-480);
				times.TSd2sTMul++;
			}
		}
		*ptdc=times;
	}	// check for last bank
}

//---------------------------------------------------------------------------------
void HandleBOR_V1190(int run, int file, int time, ITdc *ptdc)
{
	printf("\nHandleBOR_V1190...\n\n");
	if(file==0) tree->Branch("tdc","ITdc",ptdc,32000,99);
	else tree->SetBranchAddress("tdc",&ptdc);
	printf("Finished HandleBOR_V1190\n");
}

//---------------------------------------------------------------------------------
void HandleEOR_V1190(int run, int time)
{
	printf(" in V1190 EOR\n");
}
