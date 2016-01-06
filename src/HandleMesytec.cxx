//
// ROOT analyzer
//
// detector handling
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
#include <fstream>
#include <assert.h>
#include <signal.h>
#include <stdint.h>


#include "TMidasEvent.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include "TRandom3.h"
#include "TParticle.h"
#include "eloss.h"
#include <TVector3.h>
#include "runDepPar.h"

#include "HandleMesytec.h"
#include "CalibMesytec.h"
#include "../include/Globals.h"

int gMesytecnitems;

extern TEvent *IrisEvent;

extern TFile* treeFile;

extern TTree* tree;

CalibMesytec calMesy;

int usePeds = 1; //using pedestals instead of offsets for Silicon detectors AS
bool useCsI2Slope = 0; //1;
TVector3 aVector;
const int NCsI2Group = 16;
const int NCsI1GroupRing = 4;
const int Nchannels = 8*64;
const int NICChannels=16;
const int NCsIChannels=16;
const int NSd1rChannels=32;
const int NSd1sChannels=32;
const int NSd2rChannels=32;
const int NSd2sChannels=32;
const int NSurChannels=32;
const int NSusChannels=32;
const int NYdChannels=4*32;
const int NYuChannels=4*32;

char var[50];
//AS Ion Chamber
float IC[32]={0}, ICEnergy=0, ICEnergy2=0; //Dummy for IC energy
int ICChannel;  // channel with the greatest value
int ICChannel2;
float ICGain[NICChannels]={1.};
float ICPed[NICChannels]={0.};

//AS CsI
float CsI[16]={0}, CsIEnergy[16];//, CsIEnergy2; //CsI energy
int CsIChannel[16];// channel with the greatest value
// int CsIChannel2;
double CsIGain[NCsIChannels]={0.};
double CsIPed[NCsIChannels]={0.};

float CsI1[16]={0}, CsI1Energy[16];//, CsI1Energy2; //CsI energy
int CsI1Channel[16];  // channel with the greatest value
// int CsI1Channel2;
double CsI1Gain[NCsI1GroupRing][NCsIChannels]={{0.}};
double CsI1Ped[NCsIChannels]={0.};

float CsI2[16]={0}, CsI2Energy[16];//, CsI2Energy2; //CsI energy
int CsI2Channel[16];  // channel with the greatest value
// int CsI2Channel2;
double CsI2Gain[NCsI2Group][NCsIChannels]={{0.}};
double CsI2Ped[NCsIChannels]={0.};
double CsI2Slope[NCsIChannels]={0.};//YY1 ring dependence

//AS S3
float Sd2r[NSd2rChannels];
float Sd2rEnergy=0, Sd2rEnergy2= 0; //Dummy for Sd2r energy
int Sd2rChannel,  Sd2rChannel2; // channel with the greatest value
float Sd2rGain[NSd2rChannels]={1.};
float Sd2rOffset[NSd2rChannels]={0.};
float Sd2rPed[NSd1sChannels]={0.};
Bool_t S3Hit;

float Sd2s[NSd2sChannels];
float Sd2sEnergy=0,  Sd2sEnergy2= 0; //Dummy for Sd2s energy  
int Sd2sChannel, Sd2sChannel2; // channel with the greatest value
float Sd2sGain[NSd2sChannels]={1.};
float Sd2sOffset[NSd2sChannels]={0.};
float Sd2sPed[NSd1sChannels]={0.};

float Sd1r[NSd1rChannels];
float Sd1rEnergy=0, Sd1rEnergy2= 0; //Dummy for Sd1r energy
int Sd1rChannel,  Sd1rChannel2; // channel with the greatest value
float Sd1rGain[NSd1rChannels]={1.};
float Sd1rOffset[NSd1rChannels]={0.};
float Sd1rPed[NSd1sChannels]={0.};
float Sd1rGain2[NSd1rChannels]={1.}; //recalibration parameters
float Sd1rOffset2[NSd1rChannels]={0.}; //recalibration parameters

float Sd1s[NSd1sChannels];
float Sd1sEnergy=0,  Sd1sEnergy2=0; //Dummy for Sd1s energy      
int Sd1sChannel, Sd1sChannel2; // channel with the greatest value                                                                                    
float Sd1sGain[NSd1sChannels]={1.};
float Sd1sOffset[NSd1sChannels]={0.};
float Sd1sPed[NSd1sChannels]={0.};

float Sur[NSurChannels];
float SurEnergy=0, SurEnergy2= 0; //Dummy for Sur energy
int SurChannel,  SurChannel2; // channel with the greatest value
float SurGain[NSurChannels]={1.};
float SurOffset[NSurChannels]={0.};
float SurPed[NSusChannels]={0.};

float Sus[NSusChannels];
float SusEnergy=0,  SusEnergy2=0; //Dummy for Sus energy      
int SusChannel, SusChannel2; // channel with the greatest value                                                                                    
float SusGain[NSusChannels]={1.};
float SusOffset[NSusChannels]={0.};
float SusPed[NSusChannels]={0.};

float Yd[NYdChannels] ={0}; 
float YdEnergy[NYdChannels]={0.};//, YdEnergy2=0; //Dummy for Yd energy
int YdChannel[NYdChannels]={-1};//, YdChannel2; // channel with the greatest value
float YdGain[NYdChannels]={1.};
float YdOffset[NYdChannels]={0.};

float Yu[NYuChannels] ={0}; 
float YuEnergy=0, YuEnergy2=0; //Dummy for Yu energy
int YuChannel, YuChannel2; // channel with the greatest value
float YuGain[NYuChannels]={1.};
float YuOffset[NYuChannels]={0.};
//SSB
float SSBEnergy = 0;
float SSBOffset=0;
float SSBGain=0;

TGraph *g_elossSi={NULL}; // energy loss function in Si

extern FILE* ASCIIYY1;
extern FILE* ASCIICsI;
extern FILE* ASCIISd1;
extern FILE* ASCIISd2;
// extern FILE* ASCIISu;
extern FILE* ASCIIIC;
int ascii =0; // bool for writing ascii AS

uint32_t adcThresh = 0;
int ydNo;
//AS Total energies

Double_t SdETot = 0;
float YdDistance = 0.; // distance from target in mm
float Yd1r= 0., Yd2r = 0. ; // inner and outer radii in mm
float Sd1Distance = 0., Sd2Distance = 0.; //distance from target in mm
float Sdr1 = 0., Sdr2= 0.; //AS Inner and outer radii of an S3 detector (in mm).

Double_t xShift = 0;//1.97;
Double_t yShift = 0;//1.3;

float SuDistance = 200; //AS distance from target in mm
float YdTheta[128] = {0.}, SdTheta=0., phi = 0, theta1 = 0, phi1=0; //AS Dummies for theta and phi, to be used for filling histograms
TRandom3 fRandom(0);
Double_t randm; //random number between 0 and 1 for each event

uint32_t modid, oformat, vpeak, resolution, evlength, timestamp;
uint32_t channel, overflow;
 
int clearDetectors()
{
	for (int j=0; j<NICChannels; j++)	IC[j] = 0;
	for (int j=0; j<NCsIChannels; j++){
		CsI[j] = 0;
		CsI1[j] = 0;
		CsI2[j] = 0;
		CsIEnergy[j]=0;
		CsI1Energy[j]=0;
		CsI2Energy[j]=0;
		CsIChannel[j]=-1;
		CsI1Channel[j]=-1;
		CsI2Channel[j]=-1;
	}
	for (int j=0; j<NSd2rChannels; j++)	 Sd2r[j] = 0;
	for (int j=0; j<NSd2sChannels; j++)	 Sd2s[j] = 0;
	for (int j=0; j<NSd1rChannels; j++)	 Sd1r[j] = 0;
	for (int j=0; j<NSd1sChannels; j++)	 Sd1s[j] = 0;
	for (int j=0; j<NYdChannels; j++){
	 	Yd[j] = 0;
		YdEnergy[j]=0.;
		YdChannel[j]=-1;
    	YdTheta[j] = 0; 
	}
	ICEnergy=0;
	Sd2rEnergy=0;
	Sd2sEnergy=0;
	Sd1rEnergy=0;
	Sd1sEnergy=0;
	Sd2rEnergy2=0;
	Sd2sEnergy2=0;
	Sd1rEnergy2=0;
	Sd1sEnergy2=0;
	SdETot=0;
	SurEnergy=0;
	SusEnergy=0;
	SurEnergy2=0;
	SusEnergy2=0;

	ICChannel=-1;
	Sd2rChannel=-1;
	Sd2sChannel=-1;
	Sd1rChannel=-1;
	Sd1sChannel=-1;
	Sd2rChannel2=-1;
	Sd2sChannel2=-1;
	Sd1rChannel2=-1;
	Sd1sChannel2=-1;
	SurChannel=-1;
	SusChannel=-1;
	SurChannel2=-1;
	SusChannel2=-1;

	SdTheta = 0;
    phi = 0;
    theta1=0; 
    phi1=0;

 	return 0;
}

void HandleMesytec(TMidasEvent& event, void* ptr, int nitems, int MYLABEL, IDet *pdet, TString CalibFile)
{
	IDet det;
  	uint32_t *data;
  	int    i, debug = 0, debug1 = 0;
  	// int eventId = event.GetEventId();

  	data = (uint32_t *) ptr;
	    gMesytecnitems = nitems;
	 
	    // Loop over all the banks
	    if (nitems != 0) {
	       	if (debug) {
				printf("Mesytec_ Evt#:%d nitems:%d\n", event.GetSerialNumber(), nitems);     
			} 
	     
	      	for (i=0 ; i<nitems ; i++) {
	
				switch (data[i] & 0xFF000000) {
					case 0x40000000: // Header
		 		 		modid = ((data[i] & 0xFF0000) >> 16);
		 		 	 	oformat = (data[i] & 0x8000) ? 1 : 0;
		 		 	 	resolution = ((data[i] & 0x7000) >> 12);
		 		 	 	evlength = (data[i] & 0xFFF);
	
		 				if (debug) {
		  					printf("Header: id:%d of:%d res:%d el:%d\n" , modid, oformat, resolution, evlength);
	
		  				break;
					case 0xC0000000: // Trailer Event
					case 0xD0000000: // Trailer Event
					case 0xE0000000: // Trailer Event
					case 0xF0000000: // Trailer Event
		  				timestamp = data[i] & 0x3FFFFFFF;
		  				// if (debug) {
		  				//   printf("Trailer: id:%d of:%d res:%d el:%d ts:%d\n"
		  				// 	   , modid, oformat, resolution, evlength, timestamp);
		  				// }
		 
		  				break;
					case 0x04000000: // Data[I] Event
		  				channel  = ((data[i] & 0x001F0000) >> 16);
	
		  				if ((channel >= 0) && (channel<16))
		    				channel = 15 - channel; //AS Flipping channels 1-16 in shapers (due to preamp box channel flipping issue) 0->15, 1->14, ..., 15->0
		  				else if ((channel>15) && (channel <32))
		    				channel = 47 - channel; //AS Flipping channels 17-32 in shapers (due to preamp box channel flipping issue) 16->31, 17->30, ...,31->16
		  				
						overflow = (data[i] & 0x4000) ? 1 : 0;
		  				vpeak    = (data[i] & 0x1FFF);
		  				//Set a software threshold 
	
		  				//if (overflow)
		  				//vpeak = 0xFFFF;
		 
		  				if (debug1  && modid==1) printf("Evt#:%d items:%d - data[%d]: %d / 0x%x\n", event.GetSerialNumber(),Nchannels, i, data[i], data[i]); 
	
		  				if (debug1  && modid==1) printf("Data: ch:%d id:%d val:%f\n", channel, modid, (float) vpeak);
		  				// RK :  Energy Calibration 
		  				
		  				if ((modid==0) && (vpeak > adcThresh) && (vpeak<3840)){ // Why 3840? MH
		    				//AS Fill histogram
	  						IC[channel] = (float)vpeak;
	  						// ICEnergy = ((float)vpeak-ICPed[channel])*ICGain[channel];
							// if (channel==16) // change MA july3
	            			if (channel==31){
	              				SSBEnergy = float(vpeak);// *SSBGain + SSBOffset;
	  							// printf("vpeak%f",float(vpeak));
	 							// printf("ssb E is %d", SSBEnergy);
							}
		  				}
		  
		  				if (modid==1 && vpeak > adcThresh && vpeak<3840){
		    				if (channel<16){
		    					CsI1[channel] = (float)vpeak;
		    	    		}	    
			  				else if (channel>=16){
		    					CsI2[channel-16] = (float)vpeak;
		    	    		}	    
		  				}
	
		  				if (modid==2 && vpeak > adcThresh && vpeak<3840){
		 					S3Hit = 1; 	    
		 					if (!usePeds){
		    					Sd2r[channel] = Sd2rOffset[channel]+Sd2rGain[channel]*(float)vpeak;
							}
		 					else if (usePeds){
		   						Sd2r[channel] = Sd2rGain[channel]*(((float)vpeak)-Sd2rPed[channel]);
	    						} 
		  				}
		  
		  				if (modid==3 && vpeak > adcThresh && vpeak<3840){
		    				S3Hit = 1;
		    				if (!usePeds){	
		    					Sd2s[channel] = Sd2sOffset[channel]+Sd2sGain[channel]*(float)vpeak;
							}
		 					else if (usePeds){
		   						Sd2s[channel] = Sd2sGain[channel]*(((float)vpeak)-Sd2sPed[channel]);
	 							}
		   					}
	
		  				if (modid==4 && vpeak > adcThresh  && vpeak<3840){
		 					S3Hit = 1; 
	 						if (!usePeds){
		    					Sd1r[channel] = Sd1rOffset[channel]+Sd1rGain[channel]*(float)vpeak;
							}
	 						else if (usePeds){
		   						Sd1r[channel] = Sd1rGain[channel]*(((float)vpeak)-Sd1rPed[channel]);
	 							// Sd1r[channel] = Sd1rOffset2[channel]+Sd1rGain2[channel]*Sd1r[channel]; //recalibration using data for gain matching the channels
	 							}
							}	
		  				if (modid==5 && vpeak > adcThresh  && vpeak<3840){
		     				S3Hit = 1;
		     				if (!usePeds){
		    					Sd1s[channel] = Sd1sOffset[channel]+Sd1sGain[channel]*(float)vpeak;
							}
	 						else if (usePeds){
		   						Sd1s[channel] = Sd1sGain[channel]*(((float)vpeak)-Sd1sPed[channel]);
							}
						}
	 
						// Upstream S3 detector. Has to be properly implemented! 
		  				if (modid==11  && vpeak> adcThresh && vpeak<3840){
		  				  	// if (channel==0) pdet->SSB =(float)vpeak;
		  				  	//SurEnergy = ((float)vpeak-SurOffset[channel+(modid-2)*32])*SurGain[channel+(modid-2)*32];
		   					Sur[channel] = Sd1rGain[channel]*(((float)vpeak)-Sd1rPed[channel]);
		  				}
		  				    
		  				if (modid==11  && vpeak > adcThresh && vpeak<3840){
		  				  	//SusEnergy = ((float)vpeak-SusOffset[channel+(modid-2)*32])*SusGain[channel+(modid-2)*32];
		   					Sus[channel] = Sd1rGain[channel]*(((float)vpeak)-Sd1rPed[channel]);
		  				}
		  
		  				if (modid>5 && modid<10 && vpeak>adcThresh  && vpeak<3840){
		    				//YdEnergy = YdOffset[channel+(modid-6)*32]+YdGain[channel+(modid-6)*32]*(float)vpeak;
		 					//Yd[channel+(modid-6)*32]=YdEnergy;
		 					Yd[channel+(modid-6)*32]=YdOffset[channel+(modid-6)*32]+YdGain[channel+(modid-6)*32]*(float)vpeak;
		    				//printf("YdEnergy: %f, vpeak: %d, gain: %f\n",YdEnergy, vpeak, YdGain[channel+(modid-6)*32]);
		    				if (channel<16) ydNo = (modid-6)*2+1; //Yd number
		    				if (channel>15) ydNo = (modid-6)*2+2;
		  				} //modid
		  
		 
		 				// Upstream YY detector. Has to be properly implemented! 
		  				if (modid>11 && modid<16 && vpeak >adcThresh  && vpeak<3840){  
		  				  	YuEnergy = ((float)vpeak-YuOffset[channel+(modid-6)*32])*YuGain[channel+(modid-6)*32];
		  				}
	
		  				//AS angles // Is this used anywhere??
		  				if (modid==2){
		  				  	// theta = TMath::RadToDeg()*(Sdr1*(24.-channel-fRandom.Rndm())+Sdr2*(channel-fRandom.Rndm()))/24./Sd1Distance; //AS theta angle for Sd (24 - number of rings)
		  				}
		  				
		  				if (modid==3){
		  				  	phi = ((channel+fRandom.Rndm())*180./32.); //AS phi angle for Sd (32 - number of sectors)
		  				}
		  				//  if (modid>5 && modid<10){
		  				  //theta = TMath::RadToDeg()*(Yd1r*(16.-channel-0.5)+Yd2r*(channel+0.5))/16./YdDistance;
		  				  //det.TYdTheta= theta;
		  				 // cout <<  " det.TYdTheta= " <<  det.TYdTheta << endl;
		  				//}
	
		
		  				if (modid==10){
		  				  	theta1 = TMath::RadToDeg()*(Sdr1*(24.-channel-0.5)+Sdr2*(channel+0.5))/24./SuDistance; //AS theta angle for Su (24 - number of rings)
		  				}
		  				
		  				if (modid==11){
		  				  	phi1 = (channel*180./32.); //AS phi angle for Su (32 - number of sectors)
		  				}
	
		  				break;
				} // switch
	      	} // for loop
	    } // nitems != 0
	    
		//After looping over banks, fill the ascii files
	    if (modid>5 && modid<10) {// check last bank
	    	//Find the two highest elements in the array
	 		Sd2rEnergy=0; Sd2rEnergy2 =0; Sd2rChannel = -10000; Sd2rChannel2 =-10000;
	 		for (int i =0; i< NSd2rChannels;i++){
	    		if (Sd2rEnergy<Sd2r[i]){
	    		  	Sd2rEnergy2 = Sd2rEnergy;
					Sd2rChannel2 = Sd2rChannel;
	    		  	Sd2rEnergy=Sd2r[i];
	    		  	Sd2rChannel = i;
					// if (Sd2rEnergy2 > 0)
					// cout << "Path1  i: " << i << " Sd2r[i]: "<< Sd2r[i]   << " Sd2rEnergy: " << Sd2rEnergy << " Sd2rEnergy2: " << Sd2rEnergy2  << endl;
				}
				else if (Sd2rEnergy2<Sd2r[i])
				{     
					//  cout << "Path2 i: " << i << " Sd2r[i]: "<< Sd2r[i]  << " Sd2rEnergy2: " << Sd2rEnergy2  << endl;
					Sd2rEnergy2=Sd2r[i];
	 				Sd2rChannel2 = i;
				} //else if
	 		} //for
	
	    	Sd2sEnergy=0; Sd2sEnergy2 =0; Sd2sChannel = -10000; Sd2sChannel2 =-10000;
	    	for (int i =0; i< NSd2sChannels;i++){
				if (Sd2sEnergy<Sd2s[i]){
					Sd2sEnergy2 = Sd2sEnergy;
					Sd2sChannel2 = Sd2sChannel;
	      			Sd2sEnergy=Sd2s[i];
	      			Sd2sChannel = i;
				}
				else if (Sd2sEnergy2<Sd2s[i]) {
			    	Sd2sEnergy2=Sd2s[i];
					Sd2sChannel2 = i;
				}
	    	} //for
	
	  		if (ascii) 
	    		fprintf(ASCIISd2," %d  %d %d %d %d %d %d %d %d \n",
				event.GetSerialNumber(), Sd2rChannel+64, (int)Sd2rEnergy,  Sd2rChannel2+64, (int)Sd2rEnergy2,  
				Sd2sChannel+96, (int)Sd2sEnergy,  Sd2sChannel2+96, (int)Sd2sEnergy2);
		
	 		Sd1rEnergy=0; Sd1rEnergy2 =0; Sd1rChannel = -10000; Sd1rChannel2 =-10000;
	    	for (int i =0; i< NSd1rChannels;i++) {
				if (Sd1rEnergy<Sd1r[i]){
	      			Sd1rEnergy2 = Sd1rEnergy;
					Sd1rChannel2 = Sd1rChannel;
	      			Sd1rEnergy=Sd1r[i];
	      			Sd1rChannel = i;
				}
				else if (Sd1rEnergy2<Sd1r[i]){
					Sd1rEnergy2=Sd1r[i];
					Sd1rChannel2 = i;
				}
	    	}//for

	 		Sd1sEnergy=0; Sd1sEnergy2 =0; Sd1sChannel = -10000; Sd1sChannel2 =-10000;
	    	for (int i =0; i< NSd1sChannels;i++){
				if (Sd1sEnergy<Sd1s[i]){
	      			Sd1sEnergy2 = Sd1sEnergy;
					Sd1sChannel2 = Sd1sChannel;
	      			Sd1sEnergy=Sd1s[i];
	      			Sd1sChannel = i;
				}
				else if (Sd1sEnergy2<Sd1s[i]) {     
					Sd1sEnergy2=Sd1s[i];
					Sd1sChannel2 = i;
				}
	    	} //for
	
	  		if (ascii) 
	   			fprintf(ASCIISd1," %d  %d %d %d %d %d %d %d %d \n",
				event.GetSerialNumber(), Sd1rChannel+64, (int)Sd1rEnergy,  Sd1rChannel2+64, (int)Sd1rEnergy2,  
				Sd1sChannel+96, (int)Sd1sEnergy,  Sd1sChannel2+96, (int)Sd1sEnergy2);

	// 		SurEnergy=0; SurEnergy2 =0; SurChannel = -10000; SurChannel2 =-10000;
	//    	for (int i =0; i< NSurChannels;i++) {
	//			if (SurEnergy<Sur[i]){
	//      			SurEnergy2 = SurEnergy;
	//				SurChannel2 = SurChannel;
	//      			SurEnergy=Sur[i];
	//      			SurChannel = i;
	//			}
	//			else if (SurEnergy2<Sur[i]){
	//				SurEnergy2=Sur[i];
	//				SurChannel2 = i;
	//			}
	//    	}//for

	// 		SusEnergy=0; SusEnergy2 =0; SusChannel = -10000; SusChannel2 =-10000;
	//    	for (int i =0; i< NSusChannels;i++){
	//			if (SusEnergy<Sus[i]){
	//      			SusEnergy2 = SusEnergy;
	//				SusChannel2 = SusChannel;
	//      			SusEnergy=Sus[i];
	//      			SusChannel = i;
	//			}
	//			else if (SusEnergy2<Sus[i]) {     
	//				SusEnergy2=Sus[i];
	//				SusChannel2 = i;
	//			}
	//    	} //for
	
	//  		if (ascii) 
	//   			fprintf(ASCIISu," %d  %d %d %d %d %d %d %d %d \n",
	//			event.GetSerialNumber(), SurChannel+64, (int)SurEnergy,  SurChannel2+64, (int)SurEnergy2,  
	//			SusChannel+96, (int)SusEnergy,  SusChannel2+96, (int)SusEnergy2);
	//
		  	//root tree values
		    det.Clear();  
			
			det.SSB = SSBEnergy;
			det.TSd2rEnergy = Sd2rEnergy;
			det.TSd2sEnergy = Sd2sEnergy;
			det.TSd1rEnergy = Sd1rEnergy;
			det.TSd1sEnergy = Sd1sEnergy;
			det.TSd2rChannel = Sd2rChannel;
			det.TSd2sChannel = Sd2sChannel;
			det.TSd1rChannel = Sd1rChannel;
			det.TSd1sChannel = Sd1sChannel;
			//det.TSurEnergy = SurEnergy;
			//det.TSusEnergy = SusEnergy;
			//det.TSurChannel = SurChannel;
			//det.TSusChannel = SusChannel;

			//Make sure to only accept neighboring channels for energy deposition division
			if(abs(det.TSd2rChannel-det.TSd2rChannel2)==1)	det.TSd2rEnergy2 = Sd2rEnergy2;
			if(abs(det.TSd2sChannel-det.TSd2sChannel2)==1)	det.TSd2sEnergy2 = Sd2sEnergy2;
			if(abs(det.TSd1rChannel-det.TSd1rChannel2)==1)	det.TSd1rEnergy2 = Sd1rEnergy2;
			if(abs(det.TSd1sChannel-det.TSd1sChannel2)==1)	det.TSd1sEnergy2 = Sd1sEnergy2;
			//if(abs(det.TSurChannel-det.TSurChannel2)==1)	det.TSurEnergy2 = SurEnergy2;
			//if(abs(det.TSusChannel-det.TSusChannel2)==1)	det.TSusEnergy2 = SusEnergy2;
		
			det.TSd2rChannel2 = Sd2rChannel2;
			det.TSd2sChannel2 = Sd2sChannel2;
			det.TSd1rChannel2 = Sd1rChannel2;
			det.TSd1sChannel2 = Sd1sChannel2;
			//det.TSurChannel2 = SurChannel2;
			//det.TSusChannel2 = SusChannel2;
			
			randm = 0.95*fRandom.Rndm(); //random number between 0 and 0.95 for each event
			det.TSdTheta = TMath::RadToDeg()*atan((Sdr1*(24.-Sd1rChannel-randm)+Sdr2*(Sd1rChannel+randm))/24./Sd1Distance); //AS theta angle for Sd (24 - number of rings)
			
			det.TSdPhi = 180.-360.*Sd1sChannel/32.;
		 	
			//Finds the channel with the highest peak
		    // YdEnergy=0; YdChannel = -1; YdEnergy2=0; YdChannel2 =-10000;
		    // for (int i =0; i< NYdChannels;i++) {
			//   	if(Yd[i]>YdEnergy){
		  	// 		YdEnergy2= YdEnergy;
			//    		YdChannel2 = YdChannel;
			//   		YdEnergy=Yd[i];
			//   		YdChannel = i;
			//   	}	
		 	// 	else if (Yd[i]>YdEnergy2) {//if Yd[i] is between YdEnergy and YdEnergy2
			//     	YdEnergy2=Yd[i];
			//   		YdChannel2 = i;
			//     } //if
			// } //for
 
 			for (Int_t i=0;i<NYdChannels;i++){
      			Double_t max = 0.;
      			Int_t index = -1;
    			for (Int_t j=0;j<NYdChannels;j++){
      	    		if(Yd[j] > max){
      	        		max = Yd[j];
      	        		index = j;
      	    		}
        		}		
    
     			YdEnergy[i] = Yd[index];
        		YdChannel[i] = index;  
        		Yd[index] = 0.;
    		}

			//if(YdChannel>=0) YdEnergy =  (YdEnergy-YdOffset[YdChannel]) * YdGain[YdChannel];
			for(int i=0;i<NYdChannels;i++){	
				if (YdChannel[i]>=0){
					det.TYdEnergy[i] = YdEnergy[i];
					det.TYdChannel[i] = YdChannel[i];
			
					det.TYdNo[i] = int(YdChannel[i]/16);
					det.TYdRing[i] = YdChannel[i]%16;
					//here
					YdTheta[i] = TMath::RadToDeg()*atan((Yd1r*(16.-YdChannel[i]%16-randm)+Yd2r*(YdChannel[i]%16+randm))/16./YdDistance);
					det.TYdTheta[i]= YdTheta[i];
				}
			}
		 
		  	if (ascii) 
		    	fprintf(ASCIIYY1," %d  %d %d %d %d \n",event.GetSerialNumber(), YdChannel[0]+192, (int)YdEnergy[0],  YdChannel[1]+192, (int)YdEnergy[1]);
		
			//AS angles
			if (Sd1rEnergy>0){
				SdTheta = TMath::RadToDeg()*(Sdr1*(24.-Sd1rChannel-randm)+Sdr2*(Sd1rChannel+randm))/24./Sd1Distance; //AS theta angle for Sd (24 - number of rings)
			}
			 if (Sd1rEnergy>0){
			  	phi = (Sd1sChannel*180./32.); //AS phi angle for Sd (32 - number of sectors)
			}

			// for (Int_t i=0;i<NCsIChannels;i++){
      		// 	Double_t max = 0.;
      		// 	Int_t index = 0;
    		// 	for (Int_t j=0;j<NCsIChannels;j++){
      	    // 		if(CsI[j] > max){
      	    //     		max = CsI[j];
      	    //     		index = j;
      	    // 		}
        	// 	}		
    
     		// 	CsIEnergy[i] = CsI[index];
        	// 	CsIChannel[i] = index;  
        	// 	CsI[index] = 0.;
    		// }

			for (Int_t i=0;i<NCsIChannels;i++){
      			Double_t max = 0.;
      			Int_t index = 0;
    			for (Int_t j=0;j<NCsIChannels;j++){
      	    		if(CsI1[j] > max){
      	        		max = CsI1[j];
      	        		index = j;
      	    		}
        		}		
    
     			CsI1Energy[i] = CsI1[index];
        		CsI1Channel[i] = index;  
        		CsI1[index] = 0.;
    		}

			for (Int_t i=0;i<NCsIChannels;i++){
      			Double_t max = 0.;
      			Int_t index = 0;
    			for (Int_t j=0;j<NCsIChannels;j++){
      	    		if(CsI2[j] > max){
      	        		max = CsI2[j];
      	        		index = j;
      	    		}
        		}		
    
     			CsI2Energy[i] = CsI2[index];
        		CsI2Channel[i] = index;  
        		CsI2[index] = 0.;
    		}
			// CsI1Energy=0; CsI1Energy2 =0; CsI1Channel = -1; CsI1Channel2 =-1;
			// for (int i =0; i< 16;i++) {
			//   	// printf("CsI1 ch: %d, value %f\n", i, CsI1[i]);
			//   	if (CsI1Energy<CsI1[i]){
			//   		CsI1Energy2 = CsI1Energy;
			//   		CsI1Channel2 = CsI1Channel;
			//  		CsI1Energy=CsI1[i];
			//  		CsI1Channel = i;
			// 	}
			//   	else  if (CsI1Energy2<CsI1[i]){
			//       	CsI1Energy2=CsI1[i];
			//       	CsI1Channel2 = i;
			// 	}
			// } //for
			// 
			// CsI2Energy=0; CsI2Energy2 =0; CsI2Channel = -1; CsI2Channel2 =-1;
			// for (int i =0; i< 16;i++) {
			// 	//  printf("CsI2 ch: %d, value %f\n", i, CsI2[i]);
			// 	if (CsI2Energy<CsI2[i]){
			// 		CsI2Energy2 = CsI2Energy;
			// 		CsI2Channel2 = CsI2Channel;
			// 	    CsI2Energy=CsI2[i];
			// 	    CsI2Channel = i;
			// 	}
			// 	else  if (CsI2Energy2<CsI2[i]){
			// 	    CsI2Energy2=CsI2[i];
			// 	    CsI2Channel2 = i;
			// 	}
			// } //for

			for(int i=0; i<NCsIChannels; i++){
	    		det.TCsI1ADC[i] = CsI1Energy[i];
	    		det.TCsI2ADC[i] = CsI2Energy[i];
	
				if (CsI1Energy[i]>0 && CsI1Channel[i]>-1 && YdChannel>=0){
		      		int m = (YdChannel[0]%16)/(16/NCsI1GroupRing);
		      		CsI1Energy[i] = (CsI1Energy[i]-CsI1Ped[CsI1Channel[i]])*CsI1Gain[m][CsI1Channel[i]];   
		      		det.TCsI1Energy[i] = CsI1Energy[i]; 
		      		det.TCsI1Channel[i] = CsI1Channel[i];
		    	}
	
	 			// if (CsI1Energy[i]>0 && CsI1Channel[i]>-1)
	    		// 	CsI1Energy[i] = (CsI1Energy[i]-CsI1Ped[CsI1Channel[i]])*CsI1Gain[CsI1Channel[i]];
	
		  		if (CsI2Energy[i]>0.&& CsI2Energy[i] < 4000. && CsI2Channel[i]>-1 && YdChannel>=0){
					int m = (YdChannel[0]%16)/(16/NCsI2Group);
		        	CsI2Energy[i] = (CsI2Energy[i]-CsI2Ped[CsI2Channel[i]])*CsI2Gain[m][CsI2Channel[i]];
		        	det.TCsI2Energy[i] = CsI2Energy[i];
		        	det.TCsI2Channel[i] = CsI2Channel[i];
		      	}
	
	 			// if (CsI2Energy[i]>0 && CsI2Channel[i]>-1)
	    		// 	CsI2Energy[i] = (CsI2Energy[i]-CsI2Ped[CsI2Channel[i]])*CsI2Gain[CsI2Channel[i]];
	
	  			if (CsI2Energy[i]>0.){
	       			CsIEnergy[i] = CsI2Energy[i];//0.177*CsI1Energy[i]-23.3;
	 				CsIChannel[i] = CsI2Channel[i];
	     		}
	     		else if (CsI1Energy[i]>0.){
	   				CsIEnergy[i] = CsI1Energy[i];
	   				CsIChannel[i] = CsI1Channel[i];
	     		}
	
	    		det.TCsIEnergy[i] = CsIEnergy[i]; //for filling the tree
	    		det.TCsIChannel[i] = CsIChannel[i];
	    		
				// det.TCsI1Energy[i] = CsI1Energy[i]; //for filling the tree
	    		// det.TCsI1Channel[i] = CsI1Channel[i];
	
	  			// det.TCsI2Energy[i] = CsI2Energy[i]; //for filling the tree
	    		// det.TCsI2Channel[i] = CsI2Channel[i];
	
	 			if (useCsI2Slope && det.TYdRing[0] >-1)   
	     		 	{det.TCsI2Energy[i]= det.TCsI2Energy[i]+det.TYdRing[0]*CsI2Slope[det.TCsI2Channel[i]];}
			}

 			if (int(det.TCsI2Channel[0]) != det.TYdNo[0])//checking if the csi is behind YY1
   				{det.TCsI2Energy[0]=0;}

    		if (ascii)  fprintf(ASCIICsI," %d  %d %d %d %d \n",event.GetSerialNumber(), CsIChannel[0]+32, (int)CsIEnergy[0],  CsIChannel[1]+32, (int)CsIEnergy[1]);
    		
			ICEnergy=0; ICEnergy2 =0; ICChannel = -10000; ICChannel2 =-10000;
    		for (int i =0; i< NICChannels;i++) {
      			// printf("IC ch: %d, value %f\n", i, IC[i]);
      			if (ICEnergy<IC[i]){
	  				ICEnergy2 = ICEnergy;
	  				ICChannel2 = ICChannel;
          			ICEnergy=IC[i];
          			ICChannel = i;}
      			else  if (ICEnergy2<IC[i]){
          			ICEnergy2=IC[i];
          			ICChannel2 = i;
				}
    		} //for
   
    		if (ascii)  fprintf(ASCIIIC," %d  %d %d %d %d \n",event.GetSerialNumber(), ICChannel+32, (int)ICEnergy,  ICChannel2+32, (int)ICEnergy2);

    		//Use calibration values here
   
			det.TICEnergy = ICEnergy; //for filling the tree
			det.TICChannel = ICChannel;
			det.TIC00 = IC[0];
			det.TIC01 = IC[1];
			det.TIC02 = IC[2];	
			det.TIC03 = IC[3];
			det.TIC04 = IC[4];
			det.TIC05 = IC[5];
			det.TIC06 = IC[6];
			det.TIC07 = IC[7];
			det.TIC08 = IC[8];
			det.TIC09 = IC[9];
			det.TIC10 = IC[10];
			det.TIC11 = IC[11];
			det.TIC12 = IC[12];
			det.TIC13 = IC[13];
			det.TIC14 = IC[14];
			det.TIC15 = IC[15];
  			//  det.pileUp = IC[16];
  			//if (det.pileUp < 1500) det.pileUp = 0;
    		//printf("modid: %d , ICEnergy: %f \n",modid, ICEnergy);
  
  			//calculate the 
  			//det.TSd1rEnergy2 = 0; // for 40 Ar runs

  			det.TSd2rEnergyCal = erem(det.TSd1rEnergy,g_elossSi);

  			//Effect of shifting the beam by one mm
  			aVector.SetXYZ(0,0,Sd1Distance);
   			aVector.SetTheta(TMath::DegToRad()*det.TSdTheta);
   			aVector.SetPhi(TMath::DegToRad()*det.TSdPhi);

   			//aVector.SetX(aVector.X());
  			// aVector.SetY(aVector.Y());
   			aVector.SetX(aVector.X()+xShift);
   			aVector.SetY(aVector.Y()+yShift);

  			det.TSdTheta = aVector.Theta()*TMath::RadToDeg();
 			*pdet = det;
			} //last bank
  		}//nitems!=0 
// 	}// label
}

//---------------------------------------------------------------------------------
void HandleBOR_Mesytec(int run, int time, IDet* pdet, TString CalibFile)
{
	if(CalibFile=="0") printf("No calibration file specified!\n\n");
	calMesy.Load(CalibFile);
	calMesy.Print();

	ascii = calMesy.boolASCII;
	
	TFile f1(calMesy.fileELoss);    //energy loss graphs
	g_elossSi =  (TGraph*)f1.FindObjectAny("g_elossSi");  //18O energy loss in Si           
	f1.Close();
	char label[32]; //, sig[32];

// *************** Reading detector distances *************************************
	FILE * pFile;
	FILE * pwFile;
 	char buffer[32];
	
 	pwFile = fopen (Form("%s.log",CalibFile.Data()), "w");
	pFile=fopen(calMesy.fileGeometry,"r");

	if (pFile == NULL) {
		perror ("Error opening file!");
		fprintf(pwFile,"Error opening file!");
   	}	
	printf("Reading config file '%s'\n\n",calMesy.fileGeometry.Data());
	
	while (!feof(pFile))
	{
		if (!fgets(buffer,32,pFile)) break;
		printf("%s",buffer);
		
		char* val=strchr(buffer,'=');
		if (!val) printf("Missing = in input pFile, line: '%s'\n",buffer);
		*val=0;
		val++;
		if (*val==0) printf("Value missing for parameter %s",buffer);
		
		// parse float parameter (if any)
		float v;
		sscanf(val,"%f",&v);
		if (strcmp(buffer,"YDD")==0)	YdDistance = v;
		if (strcmp(buffer,"YDR1")==0)	Yd1r = v;
		if (strcmp(buffer,"YDR2")==0)	Yd2r = v;
		if (strcmp(buffer,"SD1D")==0)	Sd1Distance = v;
		if (strcmp(buffer,"SD2D")==0)	Sd2Distance = v;
		if (strcmp(buffer,"SDR1")==0)	Sdr1 = v;
		if (strcmp(buffer,"SDR2")==0)	Sdr2 = v;
	
	}
	fclose(pFile);

// ************************************************************************************

	//ASCII output file
    if (ascii) {
  		sprintf(label,"ASCIIYY1run%d.txt",gRunNumber);    
  		ASCIIYY1 = fopen(label,"w");
 		fprintf(ASCIIYY1,"Evt:  id1: vpeak1: id2: vpeak2 \n");

  		sprintf(label,"ASCIICsIrun%d.txt",gRunNumber); 
 		ASCIICsI = fopen(label,"w");
 		fprintf(ASCIICsI,"Evt:  id1: vpeak1: id2: vpeak2 \n");

 		sprintf(label,"ASCIIICrun%d.txt",gRunNumber); 
 		ASCIIIC = fopen(label,"w");
 		fprintf(ASCIICsI,"Evt:  id1: vpeak1: id2: vpeak2 \n");

  		sprintf(label,"ASCIISd2run%d.txt",gRunNumber); 
 		ASCIISd2 = fopen(label,"w");
 		fprintf(ASCIISd2,"Evt:  idr1: vpeakr1: idr2: vpeakr2:  ids1: vpeaks1: ids2: vpeaks2 \n");

 		sprintf(label,"ASCIISd1run%d.txt",gRunNumber);
		ASCIISd1 = fopen(label,"w");
 		fprintf(ASCIISd1,"Evt:  idr1: vpeakr1: idr2: vpeakr2:  ids1: vpeaks1: ids2: vpeaks2 \n");
	} //ascii AS

	treeFile->cd();
	// create a TTree
	tree = new TTree("Iris","iris data");
	// create one branch with all information from the YY1

	tree->Branch("det","IDet",pdet,32000,99);

	IrisEvent = new TEvent();
 	tree->Branch("IrisEvent","TEvent",&IrisEvent,32000,99);

// Temporary variables for calibration 
 	Int_t Chan=-10000;
	double a,b,c;
	int g; //for ringwise calibration of CsI

//************** Calibrate IC, not yet implemented! *********************************
	if(calMesy.fileIC==calMesy.installPath){
		printf("No calibration file for IC specified. Skipping IC calibration.\n\n");
	}
	else{
   		pFile = fopen(calMesy.fileIC, "r");

		if (pFile == NULL) {
			perror ("Error opening file");
			fprintf(pwFile,"Error opening file");
   		}  
 		else  {
			printf("Reading IC config file '%s'\n",calMesy.fileIC.Data());
			// Skip first line
  			fscanf(pFile,"%s",buffer);
  			fscanf(pFile,"%s",buffer);
 			fscanf(pFile,"%s",buffer);

			for (int i =0;i<16;i++  ){
    	   		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
				ICPed[Chan] = a;
				ICGain[Chan] = b;
				printf("ICPed %lf ICgain %lf\n",ICPed[Chan],ICGain[Chan]);
    	 	}
    	 	fclose (pFile);
			printf("\n");
 		}
	}
//*****************************************************************

//*************** Calibrate CsI1 ******************************
 	Chan=-10000;  
 	
 	pFile = fopen (calMesy.fileCsI1, "r");
   	if (pFile == NULL) {
		perror ("Error opening file");
		fprintf(pwFile,"Error opening file");
   	}  
 	else {
		printf("Reading config file '%s'\n",calMesy.fileCsI1.Data());
		// Skip first line
		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		//for (int i=0; i<1; i++){
		for (int i =0; i<16; i++){
			for (int j=0; j<NCsI1GroupRing; j++){
				fscanf(pFile,"%d%d%lf%lf",&Chan,&g,&a,&b);
				CsI1Ped[Chan-32] = a;
				CsI1Gain[g][Chan-32] = b;
              	printf("CsI1 calibration par: adc =%d\tc=%d\tpeds=%f\tgain=%f\n",Chan,g,a,b);
 			}//for
		}
		fclose (pFile);
		printf("\n");
 	}//else
// ******************************************************************
//
// ******************** Calibrate CsI2 *****************************
	if (useCsI2Slope)       // not using CsIslope
   	{
  		//sprintf(fileName,"%s/calib-files/CsI2PedGainSlopeS1147_11Li.txt",installPath);
		pFile = fopen (calMesy.fileCsI2,"r");

   		if (pFile == NULL) {
			perror ("Error opening file");
			fprintf(pwFile,"Error opening file");
   		}  

 		else  {
			printf("Reading config file '%s'\n",calMesy.fileCsI2.Data());
			// Skip first line
  			fscanf(pFile,"%s",buffer);
  			fscanf(pFile,"%s",buffer);
 			fscanf(pFile,"%s",buffer);
 			fscanf(pFile,"%s",buffer);

			for (int i =0;i<16;i++  ){
  				//Double_t c;
  				fscanf(pFile,"%d%lf%lf %lf",&Chan,&a,&b, &c);
				CsI2Ped[Chan] = a;
				CsI2Gain[0][Chan] = b;
				CsI2Slope[Chan] = c;
 				printf("a %lf b %lf c%lf \n",a,b,c);
      		}
     		fclose (pFile);
			printf("\n");
   		}
   	}
 	else {	 // not using CsIslope thingy so, it is reading this file

		pFile = fopen (calMesy.fileCsI2,"r");

		if (pFile == NULL) {
			perror ("Error opening file");
			fprintf(pwFile,"Error opening file");
   		}  

 		else  {
			printf("Reading config file '%s'\n",calMesy.fileCsI2.Data());
			// Skip first line
  			fscanf(pFile,"%s",buffer);
			fscanf(pFile,"%s",buffer);
			fscanf(pFile,"%s",buffer);
			fscanf(pFile,"%s",buffer);
			for(int i=0; i<16; i++){
				for (int j=0; j<NCsI2Group; j++){
       				fscanf(pFile,"%d%d%lf%lf",&Chan,&g,&a,&b);
					CsI2Ped[Chan-48] = a;  
					CsI2Gain[g][Chan-48] = b;  
					//printf("CsIPed %lf CsIgain %lf\n",a,b);
              		printf("CsI2 calibration par: adc =%d\tc=%d\tpeds=%f\tgain=%f\n",Chan,g,a,b);
     			}
			}
    		fclose (pFile);
			printf("\n");
  		}
	}//else

 	printf(" Reading CsI calibration parameters Done ....\n\n");
//************************************************************************

//**************** Calibrate Sd2 rings ****************************************
	Chan=-10000;

	pFile = fopen (calMesy.fileSd2r, "r");
	usePeds = 1;

   	if (pFile == NULL) {
		perror ("Error opening file");
		fprintf(pwFile,"Error opening file");
   	}  
 	else {
		printf("Reading Sd2r config file '%s'\n",calMesy.fileSd2r.Data());
		// Skip first line
  		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		for (int i =0;i<32;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if(!usePeds){
				Sd2rOffset[Chan-64] = a;
			}
       		else if (usePeds){
				Sd2rPed[Chan-64] = a;
				Sd2rGain[Chan-64] =  b;  
				printf("Sd2rPed %lf Sd2rgain %lf\n",a,b);
     		}
		}
     	fclose (pFile);
		printf("\n");
 	}
//************************************************************************

//**************** Calibrate Sd2 sectors ****************************************
	Chan=-10000;
   	pFile = fopen (calMesy.fileSd2s, "r");

 	if (pFile == NULL) {
		perror ("Error opening file");
		fprintf(pwFile,"Error opening file");
   	}  
 	else  {
		printf("Reading Sd2s config file '%s'\n",calMesy.fileSd2s.Data());
		// Skip first line
	  	fscanf(pFile,"%s",buffer);
	  	fscanf(pFile,"%s",buffer);
	 	fscanf(pFile,"%s",buffer);

		for (int i =0;i<32;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if (!usePeds){
				Sd2sOffset[Chan-96] = a;
			}
       		else if (usePeds){
				Sd2sPed[Chan-96] = a;
				Sd2sGain[Chan-96] = b;   
				printf("Sd2sPed %lf Sd2sgain %lf\n",a,b);
			}
     	}
     	fclose (pFile);
		printf("\n");
 	}

//************************************************************************


//**************** Calibrate Sd1 rings ****************************************
	Chan=-10000;
   	pFile = fopen(calMesy.fileSd1r, "r");

	if (pFile == NULL) {
		perror ("Error opening file");
		fprintf(pwFile,"Error opening file");
   	}  
 	else  {
		printf("Reading Sd1r config file '%s'\n",calMesy.fileSd1r.Data());
		// Skip first line
  		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		for (int i =0;i<24;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if (!usePeds)
				Sd1rOffset[Chan-128] = a;
       		else if (usePeds)
				Sd1rPed[Chan-128] = a;
				Sd1rGain[Chan-128] = b;
				printf("Sd1rPed %lf Sd1rgain %lf\n",Sd1rPed[Chan-128],Sd1rGain[Chan-128]);
     	}
     	fclose (pFile);
		printf("\n");
 	}
//************************************************************************

//**************** Calibrate Sd2 sectors ****************************************
	Chan=-10000;
   	pFile = fopen (calMesy.fileSd1s, "r");

	if (pFile == NULL) {
		perror ("Error opening file");
		fprintf(pwFile,"Error opening file");
   	}  
 	else  {
		printf("Reading Sd1s config file '%s'\n",calMesy.fileSd1s.Data());
  		// Skip first line
		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		for (int i =0;i<32;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if (!usePeds){
				Sd1sOffset[Chan-160] = a;
			}
       		else if (usePeds){
				Sd1sPed[Chan-160] = a;
				Sd1sGain[Chan-160] = b; 
				printf("Sd1sPed %lf Sd1sgain %lf\n",a,b);
 			//	if ((run>1655 && run < 1751 )|| (run >1916 && run <2035)){
   			//		Sd1sGain[Chan-128] = Sd1sGain[Chan-128]*13.3;//The gain was changed from #15 to #2 on the shaper // This has to go to the calibration file!!
     		//	}
			}
 		}
     	fclose (pFile);
		printf("\n");
	}
//************************************************************************

//**************** Calibrate Su rings ****************************************
	Chan=-10000;
	if(calMesy.fileSur==calMesy.installPath){
		printf("No calibration file for Sur specified. Skipping Sur calibration.\n");
	}
	else{
   		pFile = fopen(calMesy.fileSur, "r");

		if (pFile == NULL) {
			perror ("Error opening file");
			fprintf(pwFile,"Error opening file");
   		}  
 		else  {
			printf("Reading Sur config file '%s'\n",calMesy.fileSur.Data());
			// Skip first line
  			fscanf(pFile,"%s",buffer);
  			fscanf(pFile,"%s",buffer);
 			fscanf(pFile,"%s",buffer);

			for (int i =0;i<24;i++  ){
    	   		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
    	   		if (!usePeds)
					SurOffset[Chan-320] = a;
    	   		else if (usePeds)
					SurPed[Chan-320] = a;
					SurGain[Chan-320] = b;
					printf("SurPed %lf Surgain %lf\n",SurPed[Chan-320],SurGain[Chan-320]);
    	 	}
    	 	fclose (pFile);
			printf("\n");
 		}
	}
//************************************************************************

//**************** Calibrate Su sectors ****************************************
	Chan=-10000;
	if(calMesy.fileSus==calMesy.installPath){
		printf("No calibration file for Sus specified. Skipping Sus calibration.\n");
	}
	else{
   		pFile = fopen (calMesy.fileSus, "r");

		if (pFile == NULL) {
			perror ("Error opening file");
			fprintf(pwFile,"Error opening file");
   		}  
 		else  {
			printf("Reading Sus config file '%s'\n",calMesy.fileSus.Data());
  			// Skip first line
			fscanf(pFile,"%s",buffer);
  			fscanf(pFile,"%s",buffer);
 			fscanf(pFile,"%s",buffer);

			for (int i =0;i<32;i++  ){
    	   		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
    	   		if (!usePeds){
					SusOffset[Chan-352] = a;
				}
    	   		else if (usePeds){
					SusPed[Chan-352] = a;
					SusGain[Chan-352] = b; 
					printf("SusPed %lf Susgain %lf\n",a,b);
				}
 			}
    	 	fclose (pFile);
			printf("\n");
		}
	}
//************************************************************************

//**************** Calibrate Yd ****************************************
	Chan=-10000;
   	pFile = fopen (calMesy.fileYd, "r");

	if (pFile == NULL) {
  		perror ("Error opening file fro Yd");
  		fprintf(pwFile,"Error opening file for Yd");
	}
	else  {
		printf("Reading config file '%s'\n",calMesy.fileYd.Data());
 		// Skip first line
  		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);

  		for (int i =0;i< NYdChannels;i++  ){
    		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
			YdGain[Chan-192] = b;
 			YdOffset[Chan-192] = -1.*a*b;
 			printf("gain %lf ped %lf\t",b,a);
			if((i+1)%4==0) printf("\n");
    	}
    	fclose (pFile);
		printf("\n");
  	}
//************************************************************************

//**************** Calibrate Yu, not yet impemented!  ****************************************
	Chan=-10000;
	if(calMesy.fileYu==calMesy.installPath){
		printf("No calibration file for Yu specified. Skipping Yu calibration.\n");
	}
	else{
   		pFile = fopen (calMesy.fileYu, "r");

		if (pFile == NULL) {
  			perror ("Error opening file fro Yu");
  			fprintf(pwFile,"Error opening file for Yu");
		}
		else  {
			printf("Reading config file '%s'\n",calMesy.fileYu.Data());
 			// Skip first line
  			fscanf(pFile,"%s",buffer);
  			fscanf(pFile,"%s",buffer);
  			fscanf(pFile,"%s",buffer);

  			for (int i =0;i< NYuChannels;i++  ){
    			fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
				YuGain[Chan-384] = b;
 				YuOffset[Chan-384] = -1.*a*b;
 				printf("gain %lf ped %lf\t",b,a);
				if((i+1)%4==0) printf("\n");
    		}
    		fclose (pFile);
			printf("\n");
  		}
	}
//************************************************************************
}

void HandleEOR_Mesytec(int run, int time)
{
  printf(" in Mesytec EOR\n");
}
