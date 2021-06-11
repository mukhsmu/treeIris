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
#include <array>
#include <queue>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
//#include "midasServer.h"
#include "TMidasEvent.h"
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include "HandleV1190.h"
#include "Globals.h"

extern TFile* treeFile;

extern TTree* tree;

static uint32_t tRef=0, previousSN=-1;//, tRF =0;
//static uint32_t previousSN = -1;
int gV1190nitems;

const double timeSlope = 1; // ns/vchannel 
const uint32_t Nvchannels   = 64+128+128+64+64+64;
//const uint32_t choffset[] = {0, 64, 192, 320, 384, 448}; // The new order May 7/2015 !
const uint32_t choffset[] = {0, 64, 128, 192, 320, 384};

//Double_t timeRaw[512];
//Double_t timeRef[512];
//Double_t timeRF[512];
using namespace std;
array<queue<uint32_t>,512> timeRaw;
map<string,vector<array<double,2>>> calibrations;

uint32_t geo, evtcnt=0, tdc=0, bunchid, evtid, vchannel;
uint32_t modchannel, measure, wordcnt, evtidt, errflag, trailer;

//'Private' utility stuff
//------------------------------
struct ParsedLine {
  string option;
  string argument;
};
ParsedLine ParseLine(string &line);
map<string,vector<array<double,2>>> ParseConfigurationFile(string fileName);
vector<array<double,2>> ReadCalibrationFile(string fileName);
void CheckSize(string det,vector<array<double,2>> &cal);
//------------------------------

void HandleV1190(TMidasEvent& event, void* ptr, int nitems, int bank, ITdc* ptdc)
{
	int    i, debug = 0, debug1 = 0; 
	uint32_t *data = (uint32_t *) ptr;

	// published var for HandlePHYSICS
	gV1190nitems = nitems;
	// Loop over all the banks
	if (nitems != 0) {
   		if (debug) printf("V1190_ Evt#:%d nitems:%d bk#%d\n", event.GetSerialNumber(), nitems, bank); 
		for (i=0 ; i<nitems ; i++) {
			if (debug) 	printf("items:%d - data[%d]: %d / 0x%x\n", Nvchannels, i, data[i], data[i]);
      uint16_t globCh;
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
          //if ((globCh==3) && measure < 19500){ //TODO:Second condition??
   	    	//	tRF = measure; //RF signal
          //  if(debug) printf("\n    RF = %d\n",measure);
          //}
          //timeRaw[modchannel+choffset[bank]] = (Int_t)measure;
   	      //timeRaw[globCh] = (Int_t)measure;
          timeRaw.at(globCh).push(measure);
   					//timeRef[modchannel+choffset[bank]] = ( ((Double_t)measure - (Double_t)tRef)*timeSlope);
          //timeRef[globCh] = (((Double_t)measure - (Double_t)tRef)*timeSlope);
   					//timeRF[modchannel+choffset[bank]] = ( ((Double_t)measure - (Double_t)tRF)*timeSlope);
          //timeRF[globCh] = ( ((Double_t)measure - (Double_t)tRF)*timeSlope);
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
	
	if(bank==5){ 	// check for last bank
	  ITdc times;
    times.Clear(); //Seems unnecessary..?
    //For testing.
    //for(i=0; i<512; i++){
    //  queue<uint32_t> &ti = timeRaw.at(i);
    //  while(!ti.empty()){
    //    times.TChannel.push_back(i);
    //    times.TTDC.push_back(ti.front());
    //    ti.pop();
    //  }
    //}    

    if(tRef > 0){
    bool hasCal = false;
    hasCal = calibrations["IC"].size() == 64 ? true : false;
		for(i=0; i<64; i++){
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        times.TICTDC.push_back(tMeas);
        times.TICTChannel.push_back(i);
        times.TICTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["IC"].at(i);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TICTime.push_back(time);
        }
      }
      /*
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
      */
		}
    hasCal = calibrations["SD2S"].size() == 32 ? true : false;
	  for(i=64; i<96; i++){
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        times.TSd2sTDC.push_back(tMeas);
        times.TSd2sTChannel.push_back(95-i);
        times.TSd2sTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["SD2S"].at(95-i);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TSd2sTime.push_back(time);
        }
      }
      /*      
			if(timeRaw[i]>0.){
        if(debug) printf("Sd2s signal: channel = %d, timeRef[i] = %lf\n",i-64,timeRef[i]);
			//	times.TSd2sTDC.push_back(timeRaw[i]);
				times.TSd2sTime.push_back(timeRaw[i]);
				times.TSd2sTimeRF.push_back(timeRF[i]);
				times.TSd2sTChannel.push_back(95-i);
				times.TSd2sTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
      */
    }
    hasCal = calibrations["SD2R"].size() == 24 ? true : false;
		for(i=96; i<128; i++){ //104-127
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        if(i<104) continue;
        times.TSd2rTDC.push_back(tMeas);
        times.TSd2rTChannel.push_back(127-i);
        times.TSd2rTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["SD2R"].at(127-i);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TSd2rTime.push_back(time);
        }
      }
      /*
			if(timeRaw[i]>0.){
        if(debug) printf("Sd2r signal: channel = %d, timeRef[i] = %lf\n",i-86,timeRef[i]);
			//	times.TSd2rTDC.push_back(timeRaw[i]);
				times.TSd2rTime.push_back(timeRaw[i]);
				times.TSd2rTimeRF.push_back(timeRF[i]);
				times.TSd2rTChannel.push_back(127-i);
				times.TSd2rTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
      */
		}
    hasCal = calibrations["SD1S"].size() == 32 ? true : false;
		for(i=128; i<160; i++){
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        times.TSd1sTDC.push_back(tMeas);
        times.TSd1sTChannel.push_back(159-i);
        times.TSd1sTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["SD1S"].at(159-i);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TSd1sTime.push_back(time);
        }
      }
      /*
			if(timeRaw[i]>0.){
        if(debug) printf("Sd1s signal: channel = %d, timeRef[i] = %lf\n",i-192,timeRef[i]);
			//	times.TSd1sTDC.push_back(timeRaw[i]);
				times.TSd1sTime.push_back(timeRaw[i]);
				times.TSd1sTimeRF.push_back(timeRF[i]);
				times.TSd1sTChannel.push_back(159-i);
				times.TSd1sTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
      */
		}
    hasCal = calibrations["SD1R"].size() == 24 ? true : false;
		for(i=160; i<192; i++){ //168-191
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        if(i<168) continue;
        times.TSd1rTDC.push_back(tMeas);
        times.TSd1rTChannel.push_back(191-i);
        times.TSd1rTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["SD1R"].at(191-i);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TSd1rTime.push_back(time);
        }
      }
      /*
			if(timeRaw[i]>0.){
        if(debug) printf("Sd1r signal: channel = %d, timeRef[i] = %lf\n",i-224,timeRef[i]);
			//	times.TSd1rTDC.push_back(timeRaw[i]);
				times.TSd1rTime.push_back(timeRaw[i]);
				times.TSd1rTimeRF.push_back(timeRF[i]);
				times.TSd1rTChannel.push_back(191-i);
				times.TSd1rTMul++;
        //printf("Sd1r channel %d, traw = %lf, tRef = %u, dT = %lf, timeRef[i] = %lf\n",191-i,timeRaw[i],tRef,timeRaw[i]-tRef,timeRef[i]);
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
      */
		}
    hasCal = calibrations["YD"].size() == 128 ? true : false;
		for(i=192; i<320; i++){
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        times.TYdTDC.push_back(tMeas);
        int YdChannel = ((i-192)/16)*16+15-(i-192)%16;
        times.TYdTChannel.push_back(YdChannel);
        times.TYdTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["YD"].at(YdChannel);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TYdTime.push_back(time);
        }
      }
      /*
			if(timeRaw[i]>0.){
        if(debug) printf("Yd signal: channel = %d, timeRef[i] = %lf\n",i-320,timeRef[i]);
			  //times.TYdTDC.push_back(timeRaw[i]);
				times.TYdTime.push_back(timeRaw[i]);
				times.TYdTimeRF.push_back(timeRF[i]);
        int YdChannel = ((i-192)/16)*16+15-(i-192)%16;
				times.TYdTChannel.push_back(YdChannel);
				times.TYdTMul++;
        //printf("Yd channel %d, traw = %lf, tRef = %u, dT = %lf, timeRef[i] = %lf\n",YdChannel,timeRaw[i],tRef,timeRaw[i]-tRef,timeRef[i]);
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
      */
		}
    hasCal = calibrations["SUS"].size() == 32 ? true : false;
		for(i=320; i<352; i++){ 
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        times.TSusTDC.push_back(tMeas);
        int SusChannel = i < 336 ? 335-i : 351-i+16;
        times.TSusTChannel.push_back(SusChannel);
        times.TSusTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["SUS"].at(SusChannel);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TSusTime.push_back(time);
        }
      }
      /*
			if(timeRaw[i]>0.){
        if(debug) printf("Sus signal: channel = %d, timeRef[i] = %lf\n",i-384,timeRef[i]);
			//	times.TSusTDC.push_back(timeRaw[i]);
				times.TSusTime.push_back(timeRaw[i]);
				times.TSusTimeRF.push_back(timeRF[i]);
        int SusChannel = i < 336 ? 335-i : 351-i+16;
				times.TSusTChannel.push_back(SusChannel);
				times.TSusTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
      */
		}
    hasCal = calibrations["SUR"].size() == 24 ? true : false;
		for(i=352; i<384; i++){
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        if(i>367 && i<376) continue;
        times.TSurTDC.push_back(tMeas);
        int SurChannel = i < 368 ? 367-i : 383-i + 16;
        times.TSurTChannel.push_back(SurChannel);
        times.TSurTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["SUR"].at(SurChannel);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TSurTime.push_back(time);
        }
      }
      /*
			if(timeRaw[i]>0.){
        if(debug) printf("Sur signal: channel = %d, timeRef[i] = %lf\n",i-416,timeRef[i]);
			//	times.TSurTDC.push_back(timeRaw[i]);
				times.TSurTime.push_back(timeRaw[i]);
				times.TSurTimeRF.push_back(timeRF[i]);
        int SurChannel = i < 368 ? 367-i : 383-i + 16;
				times.TSurTChannel.push_back(SurChannel);
				times.TSurTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
      */
		}
    hasCal = calibrations["YU"].size() == 128 ? true : false;
		for(i=384; i<512; i++){
      queue<uint32_t> &ti = timeRaw.at(i);
      while(!ti.empty()){
        uint32_t tMeas = ti.front();
        ti.pop();
        times.TYuTDC.push_back(tMeas);
        int YuChannel = i < 480 ? ((i-384)/16)*16+15-(i-384)%16 : 511-i+96;
        times.TYuTChannel.push_back(YuChannel);
        times.TYuTMul++;
        if(hasCal){
          array<double,2> &coeffs =  calibrations["YU"].at(YuChannel);
          double slope = coeffs.at(0);
          double delay = coeffs.at(1);
          double time = slope*((double)tMeas-tRef-delay);
          times.TYuTime.push_back(time);
        }
      }
      /*
			if(timeRaw[i]>0.){
        if(debug) printf("Yu signal: channel = %d, timeRef[i] = %lf\n",i-448,timeRef[i]);
			//	times.TYuTDC.push_back(timeRaw[i]);
				times.TYuTime.push_back(timeRaw[i]);
				times.TYuTimeRF.push_back(timeRF[i]);
        int YuChannel = i < 480 ? ((i-384)/16)*16+15-(i-384)%16 : 511-i+96;
				times.TYuTChannel.push_back(YuChannel);
				times.TYuTMul++;
        timeRaw[i] = -1;
        timeRef[i] = -1;
        timeRF[i] = -1;
			}
      */
		}
    }//if tRef>0
		*ptdc=times;
    //tRef = 0;
	}	// check for last bank
}

//---------------------------------------------------------------------------------
void HandleBOR_V1190(int run, int file, string configFile, ITdc *ptdc)
{
	printf("\nHandleBOR_V1190...\n\n");
	if(file==0) tree->Branch("tdc","ITdc",ptdc,32000,99);
	else tree->SetBranchAddress("tdc",&ptdc);

  if(configFile != "") calibrations = ParseConfigurationFile(configFile);

	printf("Finished HandleBOR_V1190\n");
}

//---------------------------------------------------------------------------------
void HandleEOR_V1190(int run, int time)
{
	printf(" in V1190 EOR\n");
}

/**
 * Read calibration coefficients from a simple, two-column text file.
 */
vector<array<double,2>> ReadCalibrationFile(string fileName)
{
  ifstream inFile;
  string line;
  inFile.open(fileName);
  if(!inFile.is_open()){
    cout << "ReadCalibrationFile(): Error opening file \"";
    cout << fileName << "\"." << endl;
    exit(EXIT_FAILURE);
  }
  vector<array<double,2>> calibration;
  while(getline(inFile,line)){
    double slope, delay, dummy;
    if(sscanf(line.c_str(),"%lf %lf %lf",&slope,&delay,&dummy) != 2) continue;
    calibration.push_back(array<double,2>{slope,delay});
  }
  return calibration;
}

/**
 * Parse 'configuration' file. The configuration file should contain a list of
 * files from which to read time calibration coefficients. Same format as the
 * energy calibration files.
 */
map<string,vector<array<double,2>>> ParseConfigurationFile(string fileName)
{
  ifstream inFile;
  string line;
  inFile.open(fileName);
  if(!inFile.is_open()){
    cout << "ParseConfigurationFile(): Error opening file \"";
    cout << fileName << "\"." << endl;
    exit(EXIT_FAILURE);
  }
  map<string,vector<array<double,2>>> cal;
  string path = "";
  while(getline(inFile,line)){
    ParsedLine result = ParseLine(line);
    if(result.option == "PATH"){
      path = result.argument;
      continue;
    }
    string calFile = path + result.argument;
    cal[result.option] = ReadCalibrationFile(calFile);
    CheckSize(result.option,cal[result.option]);
  }
  return cal;
}

/**
 * Parse single line from the configuration file. Each line must have the format
 *   DETECTOR=/path/to/calibration/file.txt
 * where DETECTOR could be for instance YU, SD1R, SUS etc.
 */
ParsedLine ParseLine(string &line)
{
  ParsedLine result;
  result.option = line.substr(0,line.find_first_of("="));
  result.argument = line.substr(line.find_first_of("=")+1,line.length());
  return result;
}

/**
 * Check sizes of known detectors.
 */
void CheckSize(string det,vector<array<double,2>> &cal)
{
  bool ok = true;
  if(det == "YD" && cal.size() != 128) ok = false;
  if(det == "SD1R" && cal.size() != 24) ok = false;
  if(det == "SD1S" && cal.size() != 32) ok = false;
  if(det == "SD2R" && cal.size() != 24) ok = false;
  if(det == "SD2S" && cal.size() != 32) ok = false;
  if(det == "YU" && cal.size() != 128) ok = false;
  if(det == "SUR" && cal.size() != 24) ok = false;
  if(det == "SUS" && cal.size() != 32) ok = false;

  if(!ok){
    printf("Time calibration for %s has wrong size: %d\n",det.c_str(),(int)cal.size());
    exit(EXIT_FAILURE);
  }
}
