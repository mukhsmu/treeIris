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

static uint32_t tRef=0, previousSN=-1, tRF =0;
int gV1190nitems;

const double timeSlope = 1; // ns/vchannel 
const uint32_t Nvchannels   = 64+128+128+64+64+64;
//const uint32_t choffset[] = {0, 64, 192, 320, 384, 448}; // The new order May 7/2015 !
const uint32_t choffset[] = {0, 64, 128, 192, 320, 384};

Double_t timeRaw[512];
Double_t timeRef[512];
Double_t timeRF[512];

uint32_t geo, evtcnt=0, tdc=0, bunchid, evtid, vchannel;
uint32_t modchannel, measure, wordcnt, evtidt, errflag, trailer;

void HandleV1190(TMidasEvent& event, void* ptr, int nitems, int bank, ITdc* ptdc)
{
  static bool initialised = false;
  if(!initialised){ //Horrible, ugly, aargh!
    for(int i=0; i<512; i++){
      timeRaw[i] = -1;
      timeRef[i] = -1;
      timeRF[i] = -1;
    }
    initialised = true;
  }
	ITdc times;
	uint32_t *data;
	int    i, debug = 0, debug1 = 0; 
  uint32_t globCh;
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
          if(debug) printf("0x%08x: ",data[i]);
          if(debug) printf("Global TDC header. Event counter: %d, GEO: %d\n",evtcnt,geo);
					break;
				case 0x08000000:  // TDC Header
					tdc  = 0x3 & (data[i] >> 24);
					bunchid = data[i] & 0xFFF;
					evtid   = 0xFFF & (data[i] >> 12);
          if(debug) printf("  0x%08x: ",data[i]);
          if(debug) printf("TDC header. TDC: %d, evtid: %d, bunch ID: %d\n",tdc,evtid,bunchid); 
					break;  
				case 0x00000000:  // TDC measurement
          if(debug) printf("    0x%08x: ",data[i]);
					vchannel = 0x3F & (data[i] >> 19);
          //if(tdc>1) printf("Bank: %d, tdc: %d, ch: %d\n",bank,tdc,vchannel);
          if(debug) printf(" TDC data. vchannel: %d, ",vchannel);
            globCh = vchannel + (tdc/2 * 64) + choffset[bank];
			 	  	//vchannel = (vchannel/16)*16+15-vchannel%16; //Flipping vchannels 0->15; 1->14;...;15->0; 31->16; 47->32; 63;48.
  					//modchannel = (tdc/2 * 64) + vchannel;  // changed AS May 15,2013
					// modchannel = (tdc * 32) + vchannel; // original
					measure = data[i] & 0x7FFFF;
					(data[i] & 0x04000000) ? trailer = 1 : trailer = 0;
  					// Deal with delta time in reference to the STOP (tRef)
            //printf("modchannel+choffset[bank] = %d\n",modchannel+choffset[bank]);
  					//if (modchannel+choffset[bank] == 14 && event.GetSerialNumber() != previousSN) {
            if (globCh==1 && event.GetSerialNumber() != previousSN){
    					previousSN = event.GetSerialNumber();
    					// Extract tRef (ICT_@0 (1))
    					tRef = measure;
              if(debug) printf("\n    Tref = %d\n",measure);
  	  				}  
  					//if ((modchannel+choffset[bank] == 12) && measure < 19500){ //temp fix to get the 3rd RF peak
  					if ((globCh==3) && measure < 19500){ //TODO:Second condition??
   	    				tRF = measure; //RF signal
                if(debug) printf("\n    RF = %d\n",measure);
  					}
   					//timeRaw[modchannel+choffset[bank]] = (Int_t)measure;
   					timeRaw[globCh] = (Int_t)measure;
   					//timeRef[modchannel+choffset[bank]] = ( ((Double_t)measure - (Double_t)tRef)*timeSlope);
            timeRef[globCh] = (((Double_t)measure - (Double_t)tRef)*timeSlope);
   					//timeRF[modchannel+choffset[bank]] = ( ((Double_t)measure - (Double_t)tRF)*timeSlope);
   					timeRF[globCh] = ( ((Double_t)measure - (Double_t)tRF)*timeSlope);
					if (debug1) {
						printf("evt#:%d trailer:%d bk:%d tdc:%d vchannel:%d modch:%d H#:%d measure:%d/0x%x tRef:%d/0x%x dT:%d\n"
						,evtcnt, trailer, bank, tdc, vchannel, modchannel, modchannel+choffset[bank], measure, measure, tRef, tRef, measure-tRef);
					}
					break;
				case 0x18000000: // TDC Trailer
					wordcnt = data[i] & 0xFFF;
					evtidt  = 0xFFF & (data[i] >> 12);
          if(debug) printf("    0x%08x: ",data[i]);
          if(debug) printf("TDC trailer. evtid: %d, wordcnt: %d\n",evtid,wordcnt); 
					break;
        case 0x80000000: //Global trailer.
          wordcnt = data[i] & 0x001fffe0;
          wordcnt = wordcnt >> 5;
          geo = data[i] & 0x0000001f;
          if(debug) printf("  0x%08x: ",data[i]);
          if(debug) printf(", wordCnt = %d, GEO = %d\n",wordcnt,geo);
          break;
				case 0x20000000: // TDC Error
					errflag = data[i] & 0x7FFFF;
					break;
			} // switch
		} // for loop
	} // nitems!=0

	times.Clear(); //Seems unnecessary..?
	if(bank==5){ 	// check for last bank
		for(i=0; i<64; i++){
			if(timeRef[i]>0.){
        if(debug) printf("IC signal: channel = %d, timeRef[i] = %lf\n",i,timeRef[i]);
				//times.TICTDC.push_back(timeRaw[i]);
				times.TICTime.push_back(timeRef[i]);
				times.TICTimeRF.push_back(timeRF[i]);
				times.TICTChannel.push_back(i);
				times.TICTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
		}
	  for(i=64; i<96; i++){
			if(timeRef[i]>0.){
        if(debug) printf("Sd2s signal: channel = %d, timeRef[i] = %lf\n",i-64,timeRef[i]);
			//	times.TSd2sTDC.push_back(timeRaw[i]);
				times.TSd2sTime.push_back(timeRef[i]);
				times.TSd2sTimeRF.push_back(timeRF[i]);
				times.TSd2sTChannel.push_back(95-i);
				times.TSd2sTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
    }
		for(i=104; i<128; i++){
			if(timeRef[i]>0.){
        if(debug) printf("Sd2r signal: channel = %d, timeRef[i] = %lf\n",i-86,timeRef[i]);
			//	times.TSd2rTDC.push_back(timeRaw[i]);
				times.TSd2rTime.push_back(timeRef[i]);
				times.TSd2rTimeRF.push_back(timeRF[i]);
				times.TSd2rTChannel.push_back(127-i);
				times.TSd2rTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
		}
		for(i=128; i<160; i++){
			if(timeRef[i]>0.){
        if(debug) printf("Sd1s signal: channel = %d, timeRef[i] = %lf\n",i-192,timeRef[i]);
			//	times.TSd1sTDC.push_back(timeRaw[i]);
				times.TSd1sTime.push_back(timeRef[i]);
				times.TSd1sTimeRF.push_back(timeRF[i]);
				times.TSd1sTChannel.push_back(159-i);
				times.TSd1sTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
		}
		for(i=168; i<192; i++){
			if(timeRef[i]>0.){
        if(debug) printf("Sd1r signal: channel = %d, timeRef[i] = %lf\n",i-224,timeRef[i]);
			//	times.TSd1rTDC.push_back(timeRaw[i]);
				times.TSd1rTime.push_back(timeRef[i]);
				times.TSd1rTimeRF.push_back(timeRF[i]);
				times.TSd1rTChannel.push_back(191-i);
				times.TSd1rTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
		}
		for(i=192; i<320; i++){
			if(timeRef[i]>0.){
        if(debug) printf("Yd signal: channel = %d, timeRef[i] = %lf\n",i-320,timeRef[i]);
			//	times.TYdTDC.push_back(timeRaw[i]);
				times.TYdTime.push_back(timeRef[i]);
				times.TYdTimeRF.push_back(timeRF[i]);
        int YdChannel = ((i-192)/16)*16+15-(i-192)%16;
				times.TYdTChannel.push_back(YdChannel);
				times.TYdTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
		}
		for(i=320; i<351; i++){
			if(timeRef[i]>0.){
        if(debug) printf("Sus signal: channel = %d, timeRef[i] = %lf\n",i-384,timeRef[i]);
			//	times.TSusTDC.push_back(timeRaw[i]);
				times.TSusTime.push_back(timeRef[i]);
				times.TSusTimeRF.push_back(timeRF[i]);
        int SusChannel = i < 336 ? 335-i : 351-i+16;
				times.TSusTChannel.push_back(SusChannel);
				times.TSusTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
		}
		for(i=352; i<384; i++){
			if(timeRef[i]>0.){
        if(debug) printf("Sur signal: channel = %d, timeRef[i] = %lf\n",i-416,timeRef[i]);
			//	times.TSurTDC.push_back(timeRaw[i]);
				times.TSurTime.push_back(timeRef[i]);
				times.TSurTimeRF.push_back(timeRF[i]);
        int SurChannel = i < 368 ? 367-i : 383-i + 16;
				times.TSurTChannel.push_back(SurChannel);
				times.TSurTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
		}
		for(i=384; i<512; i++){
			if(timeRef[i]>0.){
        if(debug) printf("Yu signal: channel = %d, timeRef[i] = %lf\n",i-448,timeRef[i]);
			//	times.TYuTDC.push_back(timeRaw[i]);
				times.TYuTime.push_back(timeRef[i]);
				times.TYuTimeRF.push_back(timeRF[i]);
        int YuChannel = i < 480 ? ((i-384)/16)*16+15-(i-384)%16 : 511-i+96;
				times.TYuTChannel.push_back(YuChannel);
				times.TYuTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
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
