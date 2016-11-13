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
#include <vector>

#include "TMidasEvent.h"
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include "TRandom3.h"
#include <TVector3.h>
#include "TParticle.h"

#include "HandleMesytec.h"
#include "CalibMesytec.h"
#include "geometry.h"
#include "Globals.h"

int gMesytecnitems;
extern TEvent *IrisEvent;
extern TFile* treeFile;
extern TTree* tree;

CalibMesytec calMesy;
geometry geoM;

Bool_t usePeds = 0; // 1 -> using pedestals instead of offsets for Silicon detectors AS
TVector3 aVector;
const int NCsI2Group = 16;
const int NCsI1Group = 16;
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
// int CsIMul=0;
// float CsI[16]={0}, CsIEnergy[16];//, CsIEnergy2; //CsI energy
// int CsIChannel[16];// channel with the greatest value
// double CsIGain[NCsIChannels]={0.};
// double CsIPed[NCsIChannels]={0.};

int CsI1Mul=0;
int CsI1ADC[16]={0};
float CsI1[16]={0}, CsI1Energy[16];//, CsI1Energy2; //CsI energy
int CsI1Channel[16];  // channel with the greatest value
double CsI1Gain[NCsI1Group][NCsIChannels]={{1.}};
double CsI1Ped[NCsIChannels]={0.};

int CsI2Mul=0;
int CsI2ADC[16]={0};
float CsI2[16]={0}, CsI2Energy[16];//, CsI2Energy2; //CsI energy
int CsI2Channel[16];  // channel with the greatest value
double CsI2Gain[NCsI2Group][NCsIChannels]={{1.}};
double CsI2Ped[NCsIChannels]={0.};

//AS S3
int Sd1rMul=0;
int Sd1rADC[NSd1rChannels];
float Sd1r[NSd1rChannels];
float Sd1rEnergy[NSd1rChannels]={0.}; // Sd1rEnergy2= 0; //Dummy for Sd1r energy
int Sd1rChannel[NSd1rChannels]={-1}; //  Sd1rChannel2; // channel with the greatest value
float Sd1rGain[NSd1rChannels]={1.};
float Sd1rOffset[NSd1rChannels]={0.};
float Sd1rPed[NSd1sChannels]={0.};
float Sd1rGain2[NSd1rChannels]={1.}; //recalibration parameters
float Sd1rOffset2[NSd1rChannels]={0.}; //recalibration parameters

int Sd1sMul=0;
int Sd1sADC[NSd1rChannels];
float Sd1s[NSd1sChannels];
float Sd1sEnergy[NSd1sChannels]={0.}; //  Sd1sEnergy2=0; 
int Sd1sChannel[NSd1sChannels]={-1}; // Sd1sChannel2; // channel with the greatest value                                                                                    
float Sd1sGain[NSd1sChannels]={1.};
float Sd1sOffset[NSd1sChannels]={0.};
float Sd1sPed[NSd1sChannels]={0.};

int Sd2rMul=0;
int Sd2rADC[NSd1rChannels];
float Sd2r[NSd2rChannels];
float Sd2rEnergy[NSd2rChannels]={0.}; // Sd2rEnergy2= 0; //Dummy for Sd2r energy
int Sd2rChannel[NSd2rChannels]={-1}; //  Sd2rChannel2; // channel with the greatest value
float Sd2rGain[NSd2rChannels]={1.};
float Sd2rOffset[NSd2rChannels]={0.};
float Sd2rPed[NSd1sChannels]={0.};
Bool_t S3Hit;

int Sd2sMul=0;
int Sd2sADC[NSd1rChannels];
float Sd2s[NSd2sChannels];
float Sd2sEnergy[NSd2sChannels]={0.}; //  Sd2sEnergy2= 0; //Dummy for Sd2s energy  
int Sd2sChannel[NSd2sChannels]={-1}; // Sd2sChannel2; // channel with the greatest value
float Sd2sGain[NSd2sChannels]={1.};
float Sd2sOffset[NSd2sChannels]={0.};
float Sd2sPed[NSd1sChannels]={0.};

int SurMul=0;
float Sur[NSurChannels];
float SurEnergy[NSurChannels]={0.}; // SurEnergy2= 0; //Dummy for Sur energy
int SurChannel[NSurChannels]={-1}; //  SurChannel2; // channel with the greatest value
float SurGain[NSurChannels]={1.};
float SurOffset[NSurChannels]={0.};
float SurPed[NSusChannels]={0.};

int SusMul=0;
float Sus[NSusChannels];
float SusEnergy[NSusChannels]={0.}; //  SusEnergy2=0; //Dummy for Sus energy      
int SusChannel[NSusChannels]={-1}; // SusChannel2; // channel with the greatest value                                                                                    
float SusGain[NSusChannels]={1.};
float SusOffset[NSusChannels]={0.};
float SusPed[NSusChannels]={0.};

int YdMul=0;
int YdADC[NYdChannels] ={0}; 
float Yd[NYdChannels] ={0.}; 
float YdEnergy[NYdChannels]={0.};//, YdEnergy2=0; //Dummy for Yd energy
int YdChannel[NYdChannels]={-1};//, YdChannel2; // channel with the greatest value
float YdGain[NYdChannels]={1.};
float YdOffset[NYdChannels]={0.};

int YuMul=0;
float Yu[NYuChannels] ={0}; 
float YuEnergy=0, YuEnergy2=0; //Dummy for Yu energy
int YuChannel, YuChannel2; // channel with the greatest value
float YuGain[NYuChannels]={1.};
float YuOffset[NYuChannels]={0.};
//SSB
float SSBEnergy = 0;
float SSBOffset=0;
float SSBGain=0;

// Time dependent corrections
float SiTCorrFactor = 1.;
float ICTCorrFactor = 1.;

extern FILE* ASCIIYY1;
extern FILE* ASCIICsI1;
extern FILE* ASCIICsI2;
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
float YdInnerRadius= 0., YdOuterRadius = 0. ; // inner and outer radii in mm
float Sd1Distance = 0., Sd2Distance = 0.; //distance from target in mm
float SdInnerRadius = 0., SdOuterRadius= 0.; //AS Inner and outer radii of an S3 detector (in mm).

Double_t xShift = 0;//1.97;
Double_t yShift = 0;//1.3;

float SuDistance = 200; //AS distance from target in mm
float YdTheta[128] = {0.}, SdTheta=0., phi = 0, theta1 = 0, phi1=0; //AS Dummies for theta and phi, to be used for filling histograms
TRandom3 fRandom(0);
Double_t rndm; //random number between 0 and 1 for each event

uint32_t modid, oformat, vpeak, resolution, evlength, timestamp;
uint32_t channel, overflow;
 
int clearDetectors()
{
	for (int j=0; j<NICChannels; j++)	IC[j] = 0;
	// CsIMul=0;
	CsI1Mul=0;
	CsI2Mul=0;
	for (int j=0; j<NCsIChannels; j++){
	//	CsI[j] = 0;
		CsI1[j] = 0;
		CsI2[j] = 0;
		CsI1ADC[j] = 0;
		CsI2ADC[j] = 0;
	//	CsIEnergy[j]=0;
		CsI1Energy[j]=0;
		CsI2Energy[j]=0;
	//	CsIChannel[j]=-1;
		CsI1Channel[j]=-1;
		CsI2Channel[j]=-1;
	}
	Sd1rMul=0;
	for (int j=0; j<NSd1rChannels; j++){
   		Sd1r[j] = 0;
   		Sd1rADC[j] = 0;
		Sd1rEnergy[j] =0;
		Sd1rChannel[j] =-1;
	}
	Sd1sMul=0;
	for (int j=0; j<NSd1sChannels; j++){
   		Sd1s[j] = 0;
   		Sd1sADC[j] = 0;
		Sd1sEnergy[j] =0;
		Sd1sChannel[j] =-1;
	}
	Sd2rMul=0;
	for (int j=0; j<NSd2rChannels; j++){
   		Sd2r[j] = 0;
   		Sd2rADC[j] = 0;
		Sd2rEnergy[j] =0;
		Sd2rChannel[j] =-1;
	}
	Sd2sMul=0;
	for (int j=0; j<NSd2sChannels; j++){
   		Sd2s[j] = 0;
   		Sd2sADC[j] = 0;
		Sd2sEnergy[j] =0;
		Sd2sChannel[j] =-1;
	}
	YdMul=0;
	for (int j=0; j<NYdChannels; j++){
	 	Yd[j] = 0;
	 	YdADC[j] = 0;
		YdEnergy[j]=0.;
		YdChannel[j]=-1;
    	YdTheta[j] = 0; 
	}
	SurMul=0;
	for (int j=0; j<NSurChannels; j++){
   		Sur[j] = 0;
		SurEnergy[j] =0;
		SurChannel[j] =-1;
	}
	SusMul=0;
	for (int j=0; j<NSusChannels; j++){
   		Sus[j] = 0;
		SusEnergy[j] =0;
		SusChannel[j] =-1;
	}
	
	ICEnergy=0;
	SdETot=0;

	ICChannel=-1;

	SdTheta = 0;
    phi = 0;
    theta1=0; 
    phi1=0;

 	return 0;
}

void HandleMesytec(TMidasEvent& event, void* ptr, int nitems, int bank, IDet *pdet, TString CalibFile, bool gUseRaw)
{
	IDet det;
  	uint32_t *data;
  	int    i, debug = 0, debug1 = 0;

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
						//IC[channel] = (float)vpeak;
						IC[channel] = ((float)vpeak-ICPed[channel])*ICGain[channel];
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
						Sd2rADC[channel]=vpeak;		
	 					if (!usePeds){
	    					Sd2r[channel] = Sd2rOffset[channel]+Sd2rGain[channel]*(float)vpeak;
						}
	 					else if (usePeds){
	   						Sd2r[channel] = Sd2rGain[channel]*(((float)vpeak)-Sd2rPed[channel]);
							} 
	  				}
	  
	  				if (modid==3 && vpeak > adcThresh && vpeak<3840){
	    				S3Hit = 1;
						Sd2sADC[channel]=vpeak;		
	    				if (!usePeds){	
	    					Sd2s[channel] = Sd2sOffset[channel]+Sd2sGain[channel]*(float)vpeak;
						}
	 					else if (usePeds){
	   						Sd2s[channel] = Sd2sGain[channel]*(((float)vpeak)-Sd2sPed[channel]);
							}
	   					}
	
	  				if (modid==4 && vpeak > adcThresh  && vpeak<3840){
	 					S3Hit = 1; 
						Sd1rADC[channel]=vpeak;		
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
						Sd1sADC[channel]=vpeak;		
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
						YdADC[channel+(modid-6)*32]=vpeak;
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
	  				  	// theta = TMath::RadToDeg()*(SdInnerRadius*(24.-channel-fRandom.Rndm())+SdOuterRadius*(channel-fRandom.Rndm()))/24./Sd1Distance; //AS theta angle for Sd (24 - number of rings)
	  				}
	  				
	  				if (modid==3){
	  				  	phi = ((channel+fRandom.Rndm())*180./32.); //AS phi angle for Sd (32 - number of sectors)
	  				}
	  				//  if (modid>5 && modid<10){
	  				  //theta = TMath::RadToDeg()*(YdInnerRadius*(16.-channel-0.5)+YdOuterRadius*(channel+0.5))/16./YdDistance;
	  				  //det.TYdTheta= theta;
	  				 // cout <<  " det.TYdTheta= " <<  det.TYdTheta << endl;
	  				//}
	
	  				if (modid==10){
	  				  	theta1 = TMath::RadToDeg()*(SdInnerRadius*(24.-channel-0.5)+SdOuterRadius*(channel+0.5))/24./SuDistance; //AS theta angle for Su (24 - number of rings)
	  				}
	  				
	  				if (modid==11){
	  				  	phi1 = (channel*180./32.); //AS phi angle for Su (32 - number of sectors)
	  				}
	
	  				break;
			} // switch
	  	} // for loop
	} // nitems != 0
	
	//After looping over banks, fill the root tree
	
	det.Clear(); //make sure root variables are empty

	if (modid>5 && modid<10) {// check last bank

 		for (Int_t i=0;i<NSd1rChannels;i++){
			if(gUseRaw) det.TSd1rADC.push_back(Sd1rADC[i]); 
    		Double_t max = 0.;
    		Int_t index = -1;
    		for (Int_t j=0;j<NSd1rChannels;j++){
        		if(Sd1r[j] > max){
            		max = Sd1r[j];
            		index = j;
        		}
    		}		
    
    		Sd1rEnergy[i] = Sd1r[index];
			if(Sd1rEnergy[i]>0.) Sd1rMul++;
    		Sd1rChannel[i] = index;  
    		Sd1r[index] = 0.;
    	}

		det.TSd1rMul = Sd1rMul;
 		for (Int_t i=0;i<Sd1rMul;i++){
 			det.TSd1rEnergy.push_back(Sd1rEnergy[i]);
			det.TSd1rChannel.push_back(Sd1rChannel[i]);
			rndm = 0.95*fRandom.Rndm(); //random number between 0 and 0.95 for each event
			det.TSdTheta.push_back(TMath::RadToDeg()*atan((geoM.SdInnerRadius*(24.-Sd1rChannel[i]-rndm)+geoM.SdOuterRadius*(Sd1rChannel[i]+rndm))/24./geoM.Sd1Distance)); //AS theta angle for Sd (24 - number of rings)
		}
		
		for (Int_t i=0;i<NSd1sChannels;i++){
			if(gUseRaw) det.TSd1sADC.push_back(Sd1sADC[i]); 
    		Double_t max = 0.;
    		Int_t index = -1;
    		for (Int_t j=0;j<NSd1sChannels;j++){
        		if(Sd1s[j] > max){
            		max = Sd1s[j];
            		index = j;
        		}
    		}		
    
    		Sd1sEnergy[i] = Sd1s[index];
			if(Sd1sEnergy[i]>0.) Sd1sMul++;
    		Sd1sChannel[i] = index;  
    		Sd1s[index] = 0.;
    	}

		det.TSd1sMul = Sd1sMul;
		for (Int_t i=0;i<Sd1sMul;i++){
			det.TSd1sEnergy.push_back(Sd1sEnergy[i]);
			det.TSd1sChannel.push_back(Sd1sChannel[i]);
			det.TSdPhi.push_back(180.-360.*Sd1sChannel[i]/32.);
		}

  		//Effect of shifting the beam by one mm
		for(Int_t i=0; i<Sd1rMul; i++){
  			aVector.SetXYZ(0,0,geoM.Sd1Distance);
   			aVector.SetTheta(TMath::DegToRad()*det.TSdTheta.at(i));
   			if(i<Sd1sMul) aVector.SetPhi(TMath::DegToRad()*det.TSdPhi.at(i)); // if Sd1rMul>Sd1sMul, TSdPhi[i] should be zero for some events
   			aVector.SetX(aVector.X()+geoM.xShift);
   			aVector.SetY(aVector.Y()+geoM.yShift);
  			det.TSdTheta.at(i) = aVector.Theta()*TMath::RadToDeg();
		}

		if (ascii) 
			fprintf(ASCIISd1," %d  %d %d %d %d %d %d %d %d \n",
			event.GetSerialNumber(), Sd1rChannel[0]+64, (int)Sd1rEnergy[0],  Sd1rChannel[1]+64, (int)Sd1rEnergy[1],  
			Sd1sChannel[0]+96, (int)Sd1sEnergy[0],  Sd1sChannel[1]+96, (int)Sd1sEnergy[1]);

		for (Int_t i=0;i<NSd2rChannels;i++){
			if(gUseRaw) det.TSd2rADC.push_back(Sd2rADC[i]); 
    		Double_t max = 0.;
    		Int_t index = -1;
    		for (Int_t j=0;j<NSd2rChannels;j++){
        		if(Sd2r[j] > max){
            		max = Sd2r[j];
            		index = j;
        		}
    		}		
    
    		Sd2rEnergy[i] = Sd2r[index];
			if(Sd2rEnergy[i]>0.) Sd2rMul++;
    		Sd2rChannel[i] = index;  
    		Sd2r[index] = 0.;
    	}

		det.TSd2rMul = Sd2rMul;
		for (Int_t i=0;i<Sd2rMul;i++){
			det.TSd2rEnergy.push_back(Sd2rEnergy[i]);
			det.TSd2rChannel.push_back(Sd2rChannel[i]);
		}
				
		for (Int_t i=0;i<NSd2sChannels;i++){
			if(gUseRaw) det.TSd2sADC.push_back(Sd2sADC[i]); 
    		Double_t max = 0.;
    		Int_t index = -1;
    		for (Int_t j=0;j<NSd2sChannels;j++){
        		if(Sd2s[j] > max){
            		max = Sd2s[j];
            		index = j;
        		}
    		}		
    
    		Sd2sEnergy[i] = Sd2s[index];
			if(Sd2sEnergy[i]>0.) Sd2sMul++;
    		Sd2sChannel[i] = index;  
    		Sd2s[index] = 0.;
    	}

		det.TSd2sMul = Sd2sMul;
		for (Int_t i=0;i<Sd2sMul;i++){
			det.TSd2sEnergy.push_back(Sd2sEnergy[i]);
			det.TSd2sChannel.push_back(Sd2sChannel[i]);
		}
	
		if (ascii) 
			fprintf(ASCIISd2," %d  %d %d %d %d %d %d %d %d \n",
			event.GetSerialNumber(), Sd2rChannel[0]+64, (int)Sd2rEnergy[0],  Sd2rChannel[1]+64, (int)Sd2rEnergy[1],  
			Sd2sChannel[0]+96, (int)Sd2sEnergy[0],  Sd2sChannel[1]+96, (int)Sd2sEnergy[1]);

	  	//root tree values
		
		det.SSB = SSBEnergy;
		 		
		for (Int_t i=0;i<NYdChannels;i++){
			if(gUseRaw) det.TYdADC.push_back(YdADC[i]); 
    		Double_t max = 0.;
    		Int_t index = -1;
    		for (Int_t j=0;j<NYdChannels;j++){
        		if(Yd[j] > max){
            		max = Yd[j];
            		index = j;
        		}
    		}		
    
    		YdEnergy[i] = Yd[index];
			if(YdEnergy[i]>0.) YdMul++;
    		YdChannel[i] = index;  
    		Yd[index] = 0.;
    	}

		det.TYdMul = YdMul;
		for(int i=0;i<YdMul;i++){	
			det.TYdEnergy.push_back(YdEnergy[i]);
			det.TYdChannel.push_back(YdChannel[i]);
	
			det.TYdNo.push_back(int(YdChannel[i]/16));
			det.TYdRing.push_back(YdChannel[i]%16);
			//here
			rndm = fRandom.Uniform();
			YdTheta[i] = TMath::RadToDeg()*atan((geoM.YdInnerRadius*(16.-YdChannel[i]%16-rndm)+geoM.YdOuterRadius*(YdChannel[i]%16+rndm))/16./geoM.YdDistance);
			det.TYdTheta.push_back(YdTheta[i]);
		}
	 
	  	if (ascii) 
	    	fprintf(ASCIIYY1," %d  %d %d %d %d \n",event.GetSerialNumber(), YdChannel[0]+192, (int)YdEnergy[0],  YdChannel[1]+192, (int)YdEnergy[1]);
	
		for (Int_t i=0;i<NCsIChannels;i++){
			if(gUseRaw) det.TCsI1ADC.push_back(CsI1Energy[i]);
    		Double_t max = 0.;
    		Int_t index = 0;
    		for (Int_t j=0;j<NCsIChannels;j++){
        		if(CsI1[j] > max){
            		max = CsI1[j];
            		index = j;
        		}
    		}		
    
    		CsI1Energy[i] = CsI1[index];
			if(CsI1Energy[i]>0.) CsI1Mul++;
    		CsI1Channel[i] = index;  
    		CsI1[index] = 0.;
    	}

		det.TCsI1Mul = CsI1Mul;
		for(int i=0; i<CsI1Mul; i++){
	
			if (YdMul>0){
	      		int m = (YdChannel[0]%16)/(16/NCsI1Group);
	      		CsI1Energy[i] = (CsI1Energy[i]-CsI1Ped[CsI1Channel[i]])*CsI1Gain[m][CsI1Channel[i]];   
	      		det.TCsI1Energy.push_back(CsI1Energy[i]); 
	      		det.TCsI1Channel.push_back(CsI1Channel[i]);
	    	}

	 //		if (CsI1Energy[i]>0.){
	 //   		CsIEnergy[i] = CsI1Energy[i];
	 //   		CsIChannel[i] = CsI1Channel[i];
	 //		}
		}
	
		for (Int_t i=0;i<NCsIChannels;i++){
			if(gUseRaw) det.TCsI2ADC.push_back(CsI2Energy[i]);
    		Double_t max = 0.;
    		Int_t index = 0;
    		for (Int_t j=0;j<NCsIChannels;j++){
        		if(CsI2[j] > max){
            		max = CsI2[j];
            		index = j;
        		}
    		}		
    
    		CsI2Energy[i] = CsI2[index];
			if(CsI2Energy[i]>0.) CsI2Mul++;
    		CsI2Channel[i] = index;  
    		CsI2[index] = 0.;
    	}
	
		det.TCsI2Mul = CsI2Mul;
		for(int i=0; i<CsI2Mul; i++){
	  		
			if (CsI2Energy[i] < 4000. && YdMul>0){
				int m = (YdChannel[0]%16)/(16/NCsI2Group);
	        	CsI2Energy[i] = (CsI2Energy[i]-CsI2Ped[CsI2Channel[i]])*CsI2Gain[m][CsI2Channel[i]];
	        	det.TCsI2Energy.push_back(CsI2Energy[i]);
	        	det.TCsI2Channel.push_back(CsI2Channel[i]);
	      	}
	
		//	if (CsI2Energy[i]>0.){
	   	//		CsIEnergy[i] = CsI2Energy[i];//0.177*CsI1Energy[i]-23.3;
		//		CsIChannel[i] = CsI2Channel[i];
	 	//	}

			// det.TCsIEnergy.push_back(CsIEnergy[i]); //for filling the tree
			// det.TCsIChannel.push_back(CsIChannel[i]);
		
		}
//	Has to be fixed for multi-hit!!!
 		 if ((det.TCsI2Channel.size()>0&&det.TYdNo.size()>0&&int(det.TCsI2Channel.at(0)/2) != det.TYdNo.at(0))||(det.TCsI2Channel.size()>0&&det.TYdNo.size()==0))//checking if the csi is behind YY1
   		 	{det.TCsI2Energy.at(0)=0;}

    	if (ascii)  fprintf(ASCIICsI1," %d  %d %d %d %d \n",event.GetSerialNumber(), CsI1Channel[0]+32, (int)CsI1Energy[0],  CsI1Channel[1]+32, (int)CsI1Energy[1]);
    	if (ascii)  fprintf(ASCIICsI2," %d  %d %d %d %d \n",event.GetSerialNumber(), CsI2Channel[0]+32, (int)CsI2Energy[0],  CsI2Channel[1]+32, (int)CsI2Energy[1]);
    	
		ICEnergy=0; ICEnergy2 =0; ICChannel = -1; ICChannel2 =-1;
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
   
		if(ICEnergy>0.) 
		{
			det.TICEnergy = ICEnergy; //for filling the tree
			det.TICChannel = ICChannel;
		}
		for(int i=0; i<NICChannels;i++){
			if(IC[i]>0.)	det.TIC[i] = IC[i];
		}
  
 		*pdet = det;
		} //last bank
  	}//nitems!=0 
}

//---------------------------------------------------------------------------------
void HandleBOR_Mesytec(int run, int time, IDet* pdet, std::string CalibFile)
{
	if(CalibFile=="") printf("No calibration file specified!\n\n");
	calMesy.Load(CalibFile);
	calMesy.Print();

	ascii = calMesy.boolASCII;
	
	char label[32]; //, sig[32];
	geoM.ReadGeometry(calMesy.fileGeometry.data());
// ************************************************************************************

	//ASCII output file
    if (ascii) {
  		sprintf(label,"ASCIIYY1run%d.txt",gRunNumber);    
  		ASCIIYY1 = fopen(label,"w");
 		fprintf(ASCIIYY1,"Evt:  id1: vpeak1: id2: vpeak2 \n");
		
		sprintf(label,"ASCIICsI1run%d.txt",gRunNumber); 
 		ASCIICsI1 = fopen(label,"w");
 		fprintf(ASCIICsI1,"Evt:  id1: vpeak1: id2: vpeak2 \n");

  		sprintf(label,"ASCIICsI2run%d.txt",gRunNumber); 
 		ASCIICsI2 = fopen(label,"w");
 		fprintf(ASCIICsI2,"Evt:  id1: vpeak1: id2: vpeak2 \n");

 		sprintf(label,"ASCIIICrun%d.txt",gRunNumber); 
 		ASCIIIC = fopen(label,"w");
 		fprintf(ASCIIIC,"Evt:  id1: vpeak1: id2: vpeak2 \n");

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
 	Int_t Chan=-1;
	double a,b,c;
	int g; //for ringwise calibration of CsI

//************** Calibrate IC, not yet implemented! *********************************
	FILE * pFile;
	FILE * logFile;
 	char buffer[32];

	// logfile
   	logFile = fopen("treeIris.log","w");
	
   	pFile = fopen(calMesy.fileIC.data(), "r");

	if (pFile == NULL || calMesy.boolIC==false) {
		//perror ("No file");
		fprintf(logFile,"No calibration file for IC. Skipping IC calibration.\n");
		printf("No calibration file for IC. Skipping IC calibration.\n");
		for (int i =0;i<16;i++  ){
			ICPed[i] = 0.;
			ICGain[i] = 1.;
     	}
   	}  
 	else  {
		printf("Reading IC config file '%s'\n",calMesy.fileIC.data());
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
//*****************************************************************

//*************** Calibrate CsI1 ******************************
 	Chan=-1;  
 	
 	pFile = fopen (calMesy.fileCsI1.data(), "r");
   	if (pFile == NULL || calMesy.boolCsI1==false) {
		// perror ("No file");
		fprintf(logFile,"No calibration file for CsI1. Skipping CsI1 calibration.\n");
		printf("No calibration file for CsI1. Skipping CsI1 calibration.\n");
   		for (int i =0; i<16; i++){
			CsI1Ped[i] = 0.;
			for (int j=0; j<NCsI1Group; j++){
				CsI1Gain[j][i] = 1.;
 			}//for
		}
	}  
 	else {
		printf("Reading config file '%s'\n",calMesy.fileCsI1.data());
		// Skip first line
		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		//for (int i=0; i<1; i++){
		for (int i =0; i<16; i++){
			for (int j=0; j<NCsI1Group; j++){
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
	pFile = fopen (calMesy.fileCsI2.data(),"r");

	if (pFile == NULL || calMesy.boolCsI2==false) {
		fprintf(logFile,"No calibration file for CsI2. Skipping CsI2 calibration.\n");
		printf("No calibration file for CsI2. Skipping CsI2 calibration.\n");
		for (int i =0; i<16; i++){
			CsI2Ped[i] = 0.;
			for (int j=0; j<NCsI2Group; j++){
				CsI2Gain[j][i] = 1.;
 			}//for
		}
	}  

	else  {
		printf("Reading config file '%s'\n",calMesy.fileCsI2.data());
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

 	printf(" Reading CsI calibration parameters Done ....\n\n");
//************************************************************************

//**************** Calibrate Sd2 rings ****************************************
	Chan=-1;

	pFile = fopen (calMesy.fileSd2r.data(), "r");

   	if (pFile == NULL || calMesy.boolSd2r==false) {
		fprintf(logFile,"No calibration file for Sd2 rings. Skipping Sd2r calibration.\n");
		printf("No calibration file for Sd2 rings. Skipping Sd2r calibration.\n");
   		for (int i =0;i<24;i++  ){
			Sd2rPed[i] = 0.;
			Sd2rOffset[i] = 0.;
			Sd2rGain[i] = 1.;  
		}
	}  
 	else {
		printf("Reading Sd2r config file '%s'\n",calMesy.fileSd2r.data());
		// Skip first line
  		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		for (int i =0;i<24;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if(!usePeds){
				Sd2rOffset[Chan-64] = a;
				Sd2rGain[Chan-64] =  b;  
				printf("Sd2rOffset %lf Sd2rgain %lf\n",a,b);
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
	Chan=-1;
   	pFile = fopen (calMesy.fileSd2s.data(), "r");

 	if (pFile == NULL || calMesy.boolSd2s==false) {
		fprintf(logFile,"No calibration file for Sd2 sectors. Skipping Sd2s calibration.\n");
		printf("No calibration file for Sd2 sectors. Skipping Sd2s calibration.\n");
   		for (int i =0;i<32;i++  ){
			Sd2sPed[i] = 0.;
			Sd2sOffset[i] = 0.;
			Sd2sGain[i] = 1.;  
		}
	}  
 	else  {
		printf("Reading Sd2s config file '%s'\n",calMesy.fileSd2s.data());
		// Skip first line
	  	fscanf(pFile,"%s",buffer);
	  	fscanf(pFile,"%s",buffer);
	 	fscanf(pFile,"%s",buffer);

		for (int i =0;i<32;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if (!usePeds){
				Sd2sOffset[Chan-96] = a;
				Sd2sGain[Chan-96] = b;   
				printf("Sd2sOffset %lf Sd2sgain %lf\n",a,b);
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
	Chan=-1;
   	pFile = fopen(calMesy.fileSd1r.data(), "r");

	if (pFile == NULL || calMesy.boolSd1r==false) {
		fprintf(logFile,"No calibration file for Sd1 rings. Skipping Sd1r calibration.\n");
		printf("No calibration file for Sd1 rings. Skipping Sd1r calibration.\n");
   		for (int i =0;i<24;i++  ){
			Sd1rPed[i] = 0.;
			Sd1rOffset[i] = 0.;
			Sd1rGain[i] = 1.;  
		}
	}  
 	else  {
		printf("Reading Sd1r config file '%s'\n",calMesy.fileSd1r.data());
		// Skip first line
  		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		for (int i =0;i<24;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if (!usePeds){
				Sd1rOffset[Chan-128] = a;
				Sd1rGain[Chan-128] = b;
				printf("Sd1rOffset %lf Sd1rGain %lf\n",Sd1rOffset[Chan-128],Sd1rGain[Chan-128]);
			}
			else if (usePeds){
				Sd1rPed[Chan-128] = a;
				Sd1rGain[Chan-128] = b;
				printf("Sd1rPed %lf Sd1rGain %lf\n",Sd1rPed[Chan-128],Sd1rGain[Chan-128]);
			}
     	}
     	fclose (pFile);
		printf("\n");
 	}
//************************************************************************

//**************** Calibrate Sd1 sectors ****************************************
	Chan=-1;
   	pFile = fopen (calMesy.fileSd1s.data(), "r");

	if (pFile == NULL || calMesy.boolSd1s==false) {
		fprintf(logFile,"No calibration file for Sd1 sectors. Skipping Sd1s calibration.\n");
		printf("No calibration file for Sd1 sectors. Skipping Sd1s calibration.\n");
   		for (int i =0;i<32;i++  ){
			Sd1sPed[i] = 0.;
			Sd1sOffset[i] = 0.;
			Sd1sGain[i] = 1.;  
		}
	}  
 	else  {
		printf("Reading Sd1s config file '%s'\n",calMesy.fileSd1s.data());
  		// Skip first line
		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		for (int i =0;i<32;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if (!usePeds){
				Sd1sOffset[Chan-160] = a;
				Sd1sGain[Chan-160] = b; 
				printf("Sd1sOffset %lf Sd1sgain %lf\n",a,b);
			}
       		else if (usePeds){
				Sd1sPed[Chan-160] = a;
				Sd1sGain[Chan-160] = b; 
				printf("Sd1sPed %lf Sd1sgain %lf\n",a,b);
			}
 		}
     	fclose (pFile);
		printf("\n");
	}
//************************************************************************

//**************** Calibrate Su rings ****************************************
	Chan=-1;
   	pFile = fopen(calMesy.fileSur.data(), "r");

	if (pFile == NULL || calMesy.boolSur==false) {
		fprintf(logFile,"No calibration file for Su rings. Skipping Sur calibration.\n");
		printf("No calibration file for Su rings. Skipping Sur calibration.\n");
   		for (int i =0;i<24;i++  ){
			SurPed[i] = 0.;
			SurGain[i] = 1.;  
		}
	}  
 	else  {
		printf("Reading Sur config file '%s'\n",calMesy.fileSur.data());
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
//************************************************************************

//**************** Calibrate Su sectors ****************************************
	Chan=-1;
   	pFile = fopen (calMesy.fileSus.data(), "r");

	if (pFile == NULL || calMesy.boolSus==false) {
		fprintf(logFile,"No calibration file for Su sectors. Skipping Sus calibration.\n");
		printf("No calibration file for Su sectors. Skipping Sus calibration.\n");
   		for (int i =0;i<24;i++  ){
			SusPed[i] = 0.;
			SusGain[i] = 1.;  
		}
	}  
 	else  {
		printf("Reading Sus config file '%s'\n",calMesy.fileSus.data());
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
//************************************************************************

//**************** Calibrate Yd ****************************************
	Chan=-1;
   	pFile = fopen (calMesy.fileYd.data(), "r");

	if (pFile == NULL || calMesy.boolYd==false) {
		printf("No calibration file for Yd. Skipping Yd calibration.\n");
		for (int i =0;i<NYdChannels;i++  ){
			YdOffset[i] = 0.;
			YdGain[i] = 1.;  
		}
	}
	else  {
		printf("Reading config file '%s'\n",calMesy.fileYd.data());
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
	Chan=-1;
   	pFile = fopen (calMesy.fileYu.data(), "r");

	if (pFile == NULL || calMesy.boolYu==false) {
		fprintf(logFile,"No calibration file for Yu. Skipping Yu calibration.\n");
		printf("No calibration file for Yu. Skipping Yu calibration.\n");
		for (int i =0;i<NYdChannels;i++  ){
			YuOffset[i] = 0.;
			YuGain[i] = 1.;  
		}
	}
	else  {
		printf("Reading config file '%s'\n",calMesy.fileYu.data());
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
//************************************************************************

//**************** Apply time dependent correction to gains!  ****************************************
   	pFile = fopen (calMesy.fileTCorrIC.data(), "r");
	
	int run_for_corr = 0;
	
	if (pFile == NULL || calMesy.boolTCorrIC==false) {
		fprintf(logFile,"No time dependent correction for IC Energy. Skipping time dependent IC correction.\n");
		printf("No time dependent correction for IC Energy. Skipping time dependent IC correction.\n");
		ICTCorrFactor =1.;
	}
	else  {
		printf("Reading config file '%s'\n",calMesy.fileTCorrIC.data());

		while (!feof(pFile)){
    		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
			if(Chan==run){ 
				run_for_corr = Chan;
				ICTCorrFactor = b; 
			}
    	}
    	fclose (pFile);	

		if(run_for_corr==0){ 
			printf("Run %d not in list. No correction applied!\n",run);
		}
		else{
			printf("Run: %d\tIC Gain correction: %f\n\n",Chan,ICTCorrFactor);
			for(int i=0; i<NICChannels; i++){
				ICGain[i] *= 1/ICTCorrFactor;
			}
		}
  	}

	pFile = fopen (calMesy.fileTCorrSi.data(), "r");
	
	run_for_corr = 0;

	if (pFile == NULL || calMesy.boolTCorrSi==false) {
		fprintf(logFile,"No time dependent correction for Si Energy. Skipping time dependent Si correction.\n");
		printf("No time dependent correction for Si Energy. Skipping time dependent Si correction.\n");
		SiTCorrFactor =1.;
	}
	else  {
		printf("Reading config file '%s'\n",calMesy.fileTCorrSi.data());

		while (!feof(pFile)){
    		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
			if(Chan==run){
				run_for_corr = Chan;
				SiTCorrFactor = b;
			}	   
    	}
    	fclose (pFile);
		if(run_for_corr==0){ 
			printf("Run %d not in list. No correction applied!\n",run);
		}
		else{
			printf("Run: %d\tSi Gain correction: %f\n",run_for_corr,SiTCorrFactor);
			for(int i=0; i<NCsIChannels; i++){
				for(int j=0; j<NCsI1Group; j++){
					CsI1Gain[j][i] *= 1/SiTCorrFactor;
				}
				for(int j=0; j<NCsI2Group; j++){
					CsI2Gain[j][i] *= 1/SiTCorrFactor;
				}
			}
			for(int i=0; i<NSd1rChannels; i++){
				Sd1rGain[i] *= 1/SiTCorrFactor;
			}
			for(int i=0; i<NSd1sChannels; i++){
				Sd1sGain[i] *= 1/SiTCorrFactor;
			}
			for(int i=0; i<NSd2rChannels; i++){
				Sd2rGain[i] *= 1/SiTCorrFactor;
			}
			for(int i=0; i<NSd2sChannels; i++){
				Sd2sGain[i] *= 1/SiTCorrFactor;
			}
			for(int i=0; i<NSurChannels; i++){
				SurGain[i] *= 1/SiTCorrFactor;
			}
			for(int i=0; i<NSusChannels; i++){
				SusGain[i] *= 1/SiTCorrFactor;
			}
			for(int i=0; i<NYdChannels; i++){
				YdGain[i] *= 1/SiTCorrFactor;
			}
			for(int i=0; i<NYuChannels; i++){
				YuGain[i] *= 1/SiTCorrFactor;
			}
			printf("Applied correction to detector gains.\n\n");
		}
  	}

//************************************************************************


}

void HandleEOR_Mesytec(int run, int time)
{
  printf(" in Mesytec EOR\n");
}
