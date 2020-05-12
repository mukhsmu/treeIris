// ROOT analyzer
// detector handling

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

Bool_t usePeds = 1; // 1 -> using pedestals instead of offsets for Silicon detectors AS
TVector3 aVector;
const int NChannels = 512;
const int NCsI2Group = 16;
const int NCsI1Group = 16;
const int NICChannels = 16;
const int NCsIChannels = 16;
const int NSd1rChannels = 32;
const int NSd1sChannels = 32;
const int NSd2rChannels = 32;
const int NSd2sChannels = 32;
const int NSurChannels = 32;
const int NSusChannels = 32;
const int NYdChannels = 128;
const int NYuChannels = 128;
const int NTrChannels = 3;
const int NZdxChannels = 16;
const int NZdyChannels = 16;

//AS Ion Chamber
float IC[32]={0};
float ICGain[NICChannels]={1.};
float ICPed[NICChannels]={0.};

int TrADC[NTrChannels]={0};
float TrEnergy[NTrChannels]={0};
float TrGain[NTrChannels]={1.};
float TrPed[NTrChannels]={0.};
//SSB
int SSBADC = 0;
float SSBEnergy = 0;
float SSBGain=1;
float SSBPed=0;
//Scintillator
int ScADC = 0;
float ScEnergy = 0;
float ScGain=1;
float ScPed=0;

int CsI1ADC[16]={0};
float CsI1[16]={0};//, CsI1Energy2; //CsI energy
double CsI1Gain[NCsI1Group][NCsIChannels]={{1.}};
double CsI1Ped[NCsIChannels]={0.};

int CsI2ADC[16]={0};
float CsI2[16]={0};//, CsI2Energy2; //CsI energy
double CsI2Gain[NCsI2Group][NCsIChannels]={{1.}};
double CsI2Ped[NCsIChannels]={0.};

//AS S3
int Sd1rADC[NSd1rChannels];
float Sd1r[NSd1rChannels];
float Sd1rGain[NSd1rChannels]={1.};
float Sd1rOffset[NSd1rChannels]={0.};
float Sd1rPed[NSd1sChannels]={0.};
float Sd1rGain2[NSd1rChannels]={1.}; //recalibration parameters
float Sd1rOffset2[NSd1rChannels]={0.}; //recalibration parameters

int Sd1sADC[NSd1rChannels];
float Sd1s[NSd1sChannels];
float Sd1sGain[NSd1sChannels]={1.};
float Sd1sOffset[NSd1sChannels]={0.};
float Sd1sPed[NSd1sChannels]={0.};

int Sd2rADC[NSd1rChannels];
float Sd2r[NSd2rChannels];
float Sd2rGain[NSd2rChannels]={1.};
float Sd2rOffset[NSd2rChannels]={0.};
float Sd2rPed[NSd1sChannels]={0.};

int Sd2sADC[NSd1rChannels];
float Sd2s[NSd2sChannels];
float Sd2sGain[NSd2sChannels]={1.};
float Sd2sOffset[NSd2sChannels]={0.};
float Sd2sPed[NSd1sChannels]={0.};

int SurADC[NSurChannels];
float Sur[NSurChannels];
float SurGain[NSurChannels]={1.};
float SurOffset[NSurChannels]={0.};
float SurPed[NSusChannels]={0.};

int SusADC[NSusChannels];
float Sus[NSusChannels];
float SusGain[NSusChannels]={1.};
float SusOffset[NSusChannels]={0.};
float SusPed[NSusChannels]={0.};

int YdADC[NYdChannels] ={0}; 
float Yd[NYdChannels] ={0.}; 
float YdGain[NYdChannels]={1.};
float YdOffset[NYdChannels]={0.};
float YdPedestal[NYdChannels]={0.};

int YuADC[NYuChannels] ={0}; 
float Yu[NYuChannels] ={0}; 
float YuGain[NYuChannels]={1.};
float YuOffset[NYuChannels]={0.};
float YuPedestal[NYuChannels]={0.};

int ZdxADC[NZdxChannels];
float Zdx[NZdxChannels];
float ZdxGain[NZdxChannels]={1.};
float ZdxOffset[NZdxChannels]={0.};
float ZdxPed[NZdxChannels]={0.};

int ZdyADC[NZdyChannels];
float Zdy[NZdyChannels];
float ZdyGain[NZdyChannels]={1.};
float ZdyOffset[NZdyChannels]={0.};
float ZdyPed[NZdyChannels]={0.};

// Time dependent corrections
float SiTCorrFactor = 1.;
float ICTCorrFactor = 1.;
uint32_t adcThresh = 0;

int ydNo, yuNo;

//AS Total energies
float YdDistance = 0.; // distance from target in mm
float YdInnerRadius= 0., YdOuterRadius = 0. ; // inner and outer radii in mm
float Sd1Distance = 0., Sd2Distance = 0.; //distance from target in mm
float SuDistance = 0., YuDistance = 0.; //distance from target in mm
float SdInnerRadius = 0., SdOuterRadius= 0.; //AS Inner and outer radii of an S3 detector (in mm).

Double_t xShift = 0;//1.97;
Double_t yShift = 0;//1.3;

TRandom3 fRandom(0);
Double_t rndm; //random number between 0 and 1 for each event

Double_t maxE = 0.; // used for energy sorting
Int_t maxCh = -1; // used for energy sorting
Double_t phi = 0., theta = 0.; // dummy variables for angles
Double_t xpos, ypos, radius = 0.; // dummy variables for positions

uint32_t modid, oformat, vpeak, resolution, evlength, timestamp;
uint32_t channel, overflow;
 
int clearDetectors()
{
	for (int j=0; j<NICChannels; j++)	IC[j] = 0;
	for (int j=0; j<NCsIChannels; j++){
	//	CsI[j] = 0;
		CsI1[j] = 0;
		CsI2[j] = 0;
		CsI1ADC[j] = 0;
		CsI2ADC[j] = 0;
	}
	for (int j=0; j<NSd1rChannels; j++){
   		Sd1r[j] = 0;
   		Sd1rADC[j] = 0;
	}
	for (int j=0; j<NSd1sChannels; j++){
   		Sd1s[j] = 0;
   		Sd1sADC[j] = 0;
	}
	for (int j=0; j<NSd2rChannels; j++){
   		Sd2r[j] = 0;
   		Sd2rADC[j] = 0;
	}
	for (int j=0; j<NSd2sChannels; j++){
   		Sd2s[j] = 0;
   		Sd2sADC[j] = 0;
	}
	for (int j=0; j<NYdChannels; j++){
	 	Yd[j] = 0;
	 	YdADC[j] = 0;
	}
	for (int j=0; j<NSurChannels; j++){
   		Sur[j] = 0;
   		SurADC[j] = 0;
	}
	for (int j=0; j<NSusChannels; j++){
   		Sus[j] = 0;
   		SusADC[j] = 0;
	}
	for (int j=0; j<NYuChannels; j++){
	 	Yu[j] = 0;
	 	YuADC[j] = 0;
	}
	for (int j=0; j<NTrChannels; j++){
	 	TrADC[j] = 0;
		TrEnergy[j]=0.;
	}
	for (int j=0; j<NZdxChannels; j++){
	 	ZdxADC[j] = 0;
		Zdx[j]=0.;
	}
	for (int j=0; j<NZdyChannels; j++){
	 	ZdyADC[j] = 0;
		Zdy[j]=0.;
	}

	SSBEnergy=0;
	ScEnergy=0;

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
					}	
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
	
	  				if (debug1  && modid==1) printf("Evt#:%d items:%d - data[%d]: %d / 0x%x\n", event.GetSerialNumber(),NChannels, i, data[i], data[i]); 
	
	  				if (debug1  && modid==1) printf("Data: ch:%d id:%d val:%f\n", channel, modid, (float) vpeak);

					// Ionization Chamber
	  				if ((modid==0)){ //&& (vpeak > adcThresh) && (vpeak<3840)){ 
						IC[channel] = ((float)vpeak-ICPed[channel])*ICGain[channel];
	        			if (channel==18){
	          				ScADC = vpeak;
	          				ScEnergy = (float(vpeak)-ScPed) * ScGain;
						}
						if (channel==21){
	          				TrADC[0] = vpeak;
	          				TrEnergy[0] = (float(vpeak)-TrPed[0]) * TrGain[0];
						}
						if (channel==22){
	          				TrADC[1] = vpeak;
	          				TrEnergy[1] = (float(vpeak)-TrPed[1]) * TrGain[1];
						}
						if (channel==23){
	          				TrADC[2] = vpeak;
	          				TrEnergy[2] = (float(vpeak)-TrPed[2]) * TrGain[2];
						}
						if (channel==31){
	          				SSBADC = vpeak;
	          				SSBEnergy = (float(vpeak)-SSBPed) * SSBGain;
						}
	  				}
	  				
					// CsI
	  				if (modid==1){// && vpeak > adcThresh && vpeak<3840){
	    				if (channel<16){
	    					CsI1[channel] = (float)vpeak + fRandom.Uniform(-0.5,0.5);
	    					CsI1ADC[channel] = (float)vpeak;
	    	    		}	    
		  				else if (channel>=16){
	    					CsI2[channel-16] = (float)vpeak + fRandom.Uniform(-0.5,0.5);
	    					CsI2ADC[channel-16] = (float)vpeak;
	    	    		}	    
	  				}
	
					// Second downstream S3 detector, rings
	  				if (modid==2){// && vpeak > adcThresh && vpeak<3840){
						Sd2rADC[channel]=vpeak;	
						if(channel<16) ZdyADC[channel] = vpeak;
	 					if (!usePeds){
	    					Sd2r[channel] = Sd2rOffset[channel]+Sd2rGain[channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
							if(channel<16) Zdy[channel] = ZdyOffset[channel]+ZdyGain[channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
						}
	 					else if (usePeds){
	   						Sd2r[channel] = Sd2rGain[channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-Sd2rPed[channel]);
							if(channel<16) Zdy[channel] = ZdyGain[channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-ZdyPed[channel]);
							} 
	  				}
	  
					// Second downstream S3 detector, segments
	  				if (modid==3){// && vpeak > adcThresh && vpeak<3840){
						Sd2sADC[channel]=vpeak;		
						if(channel<8) ZdxADC[7-channel] = vpeak;
						if(channel>7 && channel<16) ZdxADC[channel] = vpeak;
	    				if (!usePeds){	
	    					Sd2s[channel] = Sd2sOffset[channel]+Sd2sGain[channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
							if(channel<8) Zdx[7-channel] = ZdxOffset[7-channel]+ZdxGain[7-channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
							if(channel>7 && channel<16) Zdx[channel] = ZdxOffset[channel]+ZdxGain[channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
						}
	 					else if (usePeds){
	   						Sd2s[channel] = Sd2sGain[channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-Sd2sPed[channel]);
							if(channel<8) Zdx[7-channel] = ZdxGain[7-channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-ZdxPed[7-channel]);
							if(channel>7 && channel<16) Zdx[channel] = ZdxGain[channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-ZdxPed[channel]);
							}
	   					}
	
					// First downstream S3 detector, rings
	  				if (modid==4){// && vpeak > adcThresh  && vpeak<3840){
						Sd1rADC[channel]=vpeak;		
						if (!usePeds){
	    					Sd1r[channel] = Sd1rOffset[channel]+Sd1rGain[channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
						}
						else if (usePeds){
	   						Sd1r[channel] = Sd1rGain[channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-Sd1rPed[channel]);
						}
					}	
					// First downstream S3 detector, segments
	  				if (modid==5){// && vpeak > adcThresh  && vpeak<3840){
						Sd1sADC[channel]=vpeak;		
	     				if (!usePeds){
	    					Sd1s[channel] = Sd1sOffset[channel]+Sd1sGain[channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
						}
						else if (usePeds){
	   						Sd1s[channel] = Sd1sGain[channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-Sd1sPed[channel]);
						}
					}
	
					// Upstream S3 detector, rings
	  				if (modid==10){//  && vpeak> adcThresh && vpeak<3840){
	   					SurADC[channel]=vpeak;		
						if (!usePeds){
	    					Sur[channel] = SurOffset[channel]+SurGain[channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
						}
						else if (usePeds){
	   						Sur[channel] = SurGain[channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-SurPed[channel]);
						}
					}
	  				    
					// Upstream S3 detector, segments
	  				if (modid==11){//  && vpeak > adcThresh && vpeak<3840){
	  					SusADC[channel]=vpeak;		
	     				if (!usePeds){
	    					Sus[channel] = SusOffset[channel]+SusGain[channel]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
						}
						else if (usePeds){
	   						Sus[channel] = SusGain[channel]*(((float)vpeak + fRandom.Uniform(-0.5,0.5))-SusPed[channel]);
						}
					}
	  
	 				// Downstream YY1 detector 
	  				if (modid>5 && modid<10){// && vpeak>adcThresh  && vpeak<3840){
						YdADC[channel+(modid-6)*32]=vpeak;
	 					if(!usePeds){
							Yd[channel+(modid-6)*32]=YdOffset[channel+(modid-6)*32]+YdGain[channel+(modid-6)*32]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
						}
						else{
							Yd[channel+(modid-6)*32]=YdGain[channel+(modid-6)*32]*((float)vpeak + fRandom.Uniform(-0.5,0.5) -YdPedestal[channel+(modid-6)*32]);
						}
	    				if (channel<16) ydNo = (modid-6)*2+1; //Yd number
	    				if (channel>15) ydNo = (modid-6)*2+2;
	  				}
	  
	 				// Upstream YY1 detector 
	  				if (modid>11 && modid<16){// && vpeak >adcThresh  && vpeak<3840){  
	  				  	YuADC[channel+(modid-12)*32]=vpeak;
	 					if(!usePeds){
							Yu[channel+(modid-12)*32]=YuOffset[channel+(modid-12)*32]+YuGain[channel+(modid-12)*32]*((float)vpeak + fRandom.Uniform(-0.5,0.5));
						}
						else{
							Yu[channel+(modid-12)*32]=YuGain[channel+(modid-12)*32]*((float)vpeak + fRandom.Uniform(-0.5,0.5)-YuPedestal[channel+(modid-12)*32]);
						}
	    				if (channel<16) yuNo = (modid-12)*2+1; //Yu number
	    				if (channel>15) yuNo = (modid-12)*2+2;
					}
	
	  				break;
			} // switch
	  	} // for loop
	} // nitems != 0
	
	// After looping over banks, fill the root tree
	// Detector hits are sorted by energy, then copied to the root tree
	
	det.Clear(); //make sure root variables are empty

	if(bank==5){ // check last bank

		if(gUseRaw){	// write raw ADC values to tree
 			for(Int_t i=0;i<NSd1rChannels;i++){
			   det.TSd1rADC.push_back(Sd1rADC[i]); 
			}
			for(Int_t i=0;i<NSd1sChannels;i++){
			   det.TSd1sADC.push_back(Sd1sADC[i]); 
			}
			for(Int_t i=0;i<NSd2rChannels;i++){
			   det.TSd2rADC.push_back(Sd2rADC[i]); 
			}
			for(Int_t i=0;i<NSd2sChannels;i++){
			   det.TSd2sADC.push_back(Sd2sADC[i]); 
			}
			for(Int_t i=0;i<NSurChannels;i++){
			   det.TSurADC.push_back(SurADC[i]); 
			}
			for(Int_t i=0;i<NSusChannels;i++){
			   det.TSusADC.push_back(SusADC[i]); 
			}
			for(Int_t i=0;i<NYdChannels;i++){
			   det.TYdADC.push_back(YdADC[i]); 
			}
			for(Int_t i=0;i<NYuChannels;i++){
			   det.TYuADC.push_back(YuADC[i]); 
			}
			for(Int_t i=0;i<NZdxChannels;i++){
			   det.TZdxADC.push_back(ZdxADC[i]); 
			}
			for(Int_t i=0;i<NZdyChannels;i++){
			   det.TZdyADC.push_back(ZdyADC[i]); 
			}
			for(Int_t i=0;i<NCsIChannels;i++){
			   det.TCsI1ADC.push_back(CsI1ADC[i]);
			}
			for(Int_t i=0;i<NCsIChannels;i++){
				det.TCsI2ADC.push_back(CsI2ADC[i]);
			}
			for(Int_t  i=0; i<NICChannels;i++){
				det.TICADC.push_back(IC[i]);
			}
			for(Int_t i=0; i<NTrChannels; i++){
				det.TTrADC.push_back(TrADC[i]);
			}
			det.TSSBADC = SSBADC;
			det.TScADC = ScADC;
		}

		// 1st downstream S3, ring side
 		for (Int_t i=0;i<NSd1rChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NSd1rChannels;j++){
        		if(Sd1r[j] > maxE){
            		maxE = Sd1r[j];
            		maxCh = j;
        		}
    		}		
   			if(maxE>0.){ 
				det.TSd1rMul++;
 				det.TSd1rEnergy.push_back(Sd1r[maxCh]);
				det.TSd1rChannel.push_back(maxCh);
				det.TSd1rNeighbour.push_back(-1);
				rndm = fRandom.Rndm(); //random number between 0 and 0.99 for each event
				theta = TMath::RadToDeg()*atan((geoM.SdInnerRadius*(24.-maxCh-rndm)+geoM.SdOuterRadius*(maxCh+rndm))/24./geoM.Sd1Distance);
				det.TSd1Theta.push_back(theta); //AS theta angle for Sd (24 - number of rings)
    			Sd1r[maxCh] = 0.;
			}
			else break;
    	}
				
  /*
		// correct for charge sharing. Check if two neighbouring channels have a signal, then add them up.
		for (Int_t i=0; i<det.TSd1rMul; i++){
			if(det.TSd1rMul>i+1){
				for(Int_t j=i+1; j<det.TSd1rMul; j++){
					if(TMath::Abs(det.TSd1rChannel.at(i)-det.TSd1rChannel.at(j))==1){
						det.TSd1rEnergy.at(i)+=det.TSd1rEnergy.at(j);
						det.TSd1rNeighbour.at(i)=det.TSd1rChannel.at(j);
						det.TSd1rMul--;
						det.TSd1rEnergy.erase(det.TSd1rEnergy.begin()+j);
						det.TSd1rChannel.erase(det.TSd1rChannel.begin()+j);
						det.TSd1rNeighbour.erase(det.TSd1rNeighbour.begin()+j);
						det.TSd1Theta.erase(det.TSd1Theta.begin()+j);
						break;
					}
				}
			}
			if(det.TSd1rMul>=i) break;
		}
	*/
		// 1st downstream S3, sector side
		for (Int_t i=0;i<NSd1sChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NSd1sChannels;j++){
        		if(Sd1s[j] > maxE){
            		maxE = Sd1s[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TSd1sMul++;
    			det.TSd1sEnergy.push_back(Sd1s[maxCh]);
				det.TSd1sChannel.push_back(maxCh);
				det.TSd1sNeighbour.push_back(-1);
				phi = -180.+360.*maxCh/32.;
				rndm = fRandom.Rndm(); //random number between 0 and 0.99 for each event
				det.TSd1Phi.push_back(phi+11.25*rndm);
	
				Sd1s[maxCh] = 0.;
			}
			else break;
    	}
  /*
		// correct for charge sharing. Check if two neighbouring channels have a signal, then add them up.
		for (Int_t i=0; i<det.TSd1sMul; i++){
			if(det.TSd1sMul>i+1){
				for(Int_t j=i+1; j<det.TSd1sMul; j++){
					if(TMath::Abs(det.TSd1sChannel.at(i)-det.TSd1sChannel.at(j))==1||TMath::Abs(det.TSd1sChannel.at(i)-det.TSd1sChannel.at(j))==31){
						det.TSd1sEnergy.at(i)+=det.TSd1sEnergy.at(j);
						det.TSd1sNeighbour.at(i)=det.TSd1sChannel.at(j);
						det.TSd1sMul--;
						det.TSd1sEnergy.erase(det.TSd1sEnergy.begin()+j);
						det.TSd1sChannel.erase(det.TSd1sChannel.begin()+j);
						det.TSd1sNeighbour.erase(det.TSd1sNeighbour.begin()+j);
						det.TSd1Phi.erase(det.TSd1Phi.begin()+j);
						break;
					}
				}
			}
			if(det.TSd1sMul>=i) break;
		}
  */

		// 2nd downstream S3, ring side
		for (Int_t i=0;i<NSd2rChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NSd2rChannels;j++){
        		if(Sd2r[j] > maxE){
            		maxE = Sd2r[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TSd2rMul++;
    			det.TSd2rEnergy.push_back(Sd2r[maxCh]);
				det.TSd2rChannel.push_back(maxCh);
				det.TSd2rNeighbour.push_back(-1);
				rndm = fRandom.Rndm(); //random number between 0 and 0.99 for each event
				theta = TMath::RadToDeg()*atan((geoM.SdInnerRadius*(24.-maxCh-rndm)+geoM.SdOuterRadius*(maxCh+rndm))/24./(geoM.Sd1Distance+14.8));
				det.TSd2Theta.push_back(theta); //AS theta angle for Sd (24 - number of rings)

				Sd2r[maxCh] = 0.;
			}
			else break;
    	}

  /*			
		// correct for charge sharing. Check if two neighbouring channels have a signal, then add them up.
		for (Int_t i=0; i<det.TSd2rMul; i++){
			if(det.TSd2rMul>i+1){
				for(Int_t j=i+1; j<det.TSd2rMul; j++){
					if(TMath::Abs(det.TSd2rChannel.at(i)-det.TSd2rChannel.at(j))==1){
						det.TSd2rEnergy.at(i)+=det.TSd2rEnergy.at(j);
						det.TSd2rNeighbour.at(i)=det.TSd2rChannel.at(j);
						det.TSd2rMul--;
						det.TSd2rEnergy.erase(det.TSd2rEnergy.begin()+j);
						det.TSd2rChannel.erase(det.TSd2rChannel.begin()+j);
						det.TSd2rNeighbour.erase(det.TSd2rNeighbour.begin()+j);
						det.TSd2Theta.erase(det.TSd2Theta.begin()+j);
						break;
					}
				}
			}
			if(det.TSd2rMul>=i) break;
		}
  */		

		// 2nd downstream S3, sector side
		for (Int_t i=0;i<NSd2sChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NSd2sChannels;j++){
        		if(Sd2s[j] > maxE){
            		maxE = Sd2s[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TSd2sMul++;
    			det.TSd2sEnergy.push_back(Sd2s[maxCh]);
				det.TSd2sChannel.push_back(maxCh);
				det.TSd2sNeighbour.push_back(-1);
				phi = 180.-360.*maxCh/32.;
				rndm = fRandom.Rndm(); //random number between 0 and 0.99 for each event
				det.TSd2Phi.push_back(phi-11.25*rndm);

				Sd2s[maxCh] = 0.;
			}
			else break;
    	}

  /*		
		// correct for charge sharing. Check if two neighbouring channels have a signal, then add them up.
		for (Int_t i=0; i<det.TSd2sMul; i++){
			if(det.TSd2sMul>i+1){
				for(Int_t j=i+1; j<det.TSd2sMul; j++){
					if(TMath::Abs(det.TSd2sChannel.at(i)-det.TSd2sChannel.at(j))==1||TMath::Abs(det.TSd2sChannel.at(i)-det.TSd2sChannel.at(j))==31){
						det.TSd2sEnergy.at(i)+=det.TSd2sEnergy.at(j);
						det.TSd2sNeighbour.at(i)=det.TSd2sChannel.at(j);
						det.TSd2sMul--;
						det.TSd2sEnergy.erase(det.TSd2sEnergy.begin()+j);
						det.TSd2sChannel.erase(det.TSd2sChannel.begin()+j);
						det.TSd2sNeighbour.erase(det.TSd2sNeighbour.begin()+j);
						det.TSd2Phi.erase(det.TSd2Phi.begin()+j);
						break;
					}
				}
			}
			if(det.TSd2sMul>=i) break;
		}
  */
		// Upstream S3, ring side
		for (Int_t i=0;i<NSurChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NSurChannels;j++){
        		if(Sur[j] > maxE){
            		maxE = Sur[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TSurMul++;
    			det.TSurEnergy.push_back(Sur[maxCh]);
				det.TSurChannel.push_back(maxCh);
				det.TSurNeighbour.push_back(-1);
				rndm = fRandom.Rndm(); //random number between 0 and 0.99 for each event
				theta = TMath::RadToDeg()*atan((geoM.SdInnerRadius*(maxCh+rndm)+geoM.SdOuterRadius*(24-maxCh-rndm))/24./(geoM.SuDistance)) + 180.;
				det.TSuTheta.push_back(theta); //AS theta angle for Sd (24 - number of rings)

				Sur[maxCh] = 0.;
			}
			else break;
    	}
	
  /*		
		// correct for charge sharing. Check if two neighbouring channels have a signal, then add them up.
		for (Int_t i=0; i<det.TSurMul; i++){
			if(det.TSurMul>i+1){
				for(Int_t j=i+1; j<det.TSurMul; j++){
					if(TMath::Abs(det.TSurChannel.at(i)-det.TSurChannel.at(j))==1){
						det.TSurEnergy.at(i)+=det.TSurEnergy.at(j);
						det.TSurNeighbour.at(i)=det.TSurChannel.at(j);
						det.TSurMul--;
						det.TSurEnergy.erase(det.TSurEnergy.begin()+j);
						det.TSurChannel.erase(det.TSurChannel.begin()+j);
						det.TSurNeighbour.erase(det.TSurNeighbour.begin()+j);
						det.TSuTheta.erase(det.TSuTheta.begin()+j);
						break;
					}
				}
			}
			if(det.TSurMul>=i) break;
		}
  */

		// Upstream S3, sector side
		for (Int_t i=0;i<NSusChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NSusChannels;j++){
        		if(Sus[j] > maxE){
            		maxE = Sus[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TSusMul++;
    			det.TSusEnergy.push_back(Sus[maxCh]);
				det.TSusChannel.push_back(maxCh);
				det.TSusNeighbour.push_back(-1);
				phi = 180.-360.*maxCh/32.;
				rndm = fRandom.Rndm(); //random number between 0 and 0.99 for each event
				det.TSuPhi.push_back(phi-11.25*rndm);
	
				Sus[maxCh] = 0.;
			}
			else break;
    	}
		 /*
		// correct for charge sharing. Check if two neighbouring channels have a signal, then add them up.
		for (Int_t i=0; i<det.TSusMul; i++){
			if(det.TSusMul>i+1){
				for(Int_t j=i+1; j<det.TSusMul; j++){
					if(TMath::Abs(det.TSusChannel.at(i)-det.TSusChannel.at(j))==1||TMath::Abs(det.TSusChannel.at(i)-det.TSusChannel.at(j))==31){
						det.TSusEnergy.at(i)+=det.TSusEnergy.at(j);
						det.TSusNeighbour.at(i)=det.TSusChannel.at(j);
						det.TSusMul--;
						det.TSusEnergy.erase(det.TSusEnergy.begin()+j);
						det.TSusChannel.erase(det.TSusChannel.begin()+j);
						det.TSusNeighbour.erase(det.TSusNeighbour.begin()+j);
						det.TSuPhi.erase(det.TSuPhi.begin()+j);
						break;
					}
				}
			}
			if(det.TSusMul>=i) break;
		}
  */

		// Downstream YY1
		for (Int_t i=0;i<NYdChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NYdChannels;j++){
        		if(Yd[j] > maxE){
            		maxE = Yd[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TYdMul++;
    			det.TYdEnergy.push_back(Yd[maxCh]);
				det.TYdChannel.push_back(maxCh);
				det.TYdNeighbour.push_back(-1);
				det.TYdNo.push_back(int(maxCh/16));
				det.TYdRing.push_back(maxCh%16);
				rndm = fRandom.Rndm();
				theta = TMath::RadToDeg()*atan((geoM.YdInnerRadius*(16.-maxCh%16-rndm)+geoM.YdOuterRadius*(maxCh%16+rndm))/16./geoM.YdDistance);
				det.TYdTheta.push_back(theta);

				Yd[maxCh] = 0.;
			}
			else break;
    	}

  /*
		// correct for charge sharing. Check if two neighbouring channels have a signal, then add them up.
		for (Int_t i=0; i<det.TYdMul; i++){
			if(det.TYdMul>i+1){
				for(Int_t j=i+1; j<det.TYdMul; j++){
					if(TMath::Abs(det.TYdChannel.at(i)-det.TYdChannel.at(j))==1&&det.TYdNo.at(i)-det.TYdNo.at(j)==0){
						det.TYdEnergy.at(i)+=det.TYdEnergy.at(j);
						det.TYdNeighbour.at(i)=det.TYdChannel.at(j);
						det.TYdMul--;
						det.TYdEnergy.erase(det.TYdEnergy.begin()+j);
						det.TYdChannel.erase(det.TYdChannel.begin()+j);
						det.TYdNeighbour.erase(det.TYdNeighbour.begin()+j);
						det.TYdNo.erase(det.TYdNo.begin()+j);
						det.TYdRing.erase(det.TYdRing.begin()+j);
						det.TYdTheta.erase(det.TYdTheta.begin()+j);
						break;
					}
				}
			}
			if(det.TYdMul>=i) break;
		}
  */

		// Upstream YY1
		for (Int_t i=0;i<NYuChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NYuChannels;j++){
        		if(Yu[j] > maxE){
            		maxE = Yu[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TYuMul++;
    			det.TYuEnergy.push_back(Yu[maxCh]);
				det.TYuChannel.push_back(maxCh);
				det.TYuNeighbour.push_back(-1);
				det.TYuNo.push_back(int(maxCh/16));
				det.TYuRing.push_back(maxCh%16);
				//here
				rndm = fRandom.Rndm();
				theta = TMath::RadToDeg()*atan((geoM.YdInnerRadius*(16.-maxCh%16-rndm)+geoM.YdOuterRadius*(maxCh%16+rndm))/16./geoM.YuDistance) + 180.;
				det.TYuTheta.push_back(theta);

				Yu[maxCh] = 0.;
			}
			else break;
    	}
	
  /*
		// correct for charge sharing. Check if two neighbouring channels have a signal, then add them up.
		for (Int_t i=0; i<det.TYuMul; i++){
			if(det.TYuMul>i+1){
				for(Int_t j=i+1; j<det.TYuMul; j++){
					if(TMath::Abs(det.TYuChannel.at(i)-det.TYuChannel.at(j))==1&&det.TYuNo.at(i)-det.TYuNo.at(j)==0.){
						det.TYuEnergy.at(i)+=det.TYuEnergy.at(j);
						det.TYuNeighbour.at(i)=det.TYuChannel.at(j);
						det.TYuMul--;
						det.TYuEnergy.erase(det.TYuEnergy.begin()+j);
						det.TYuChannel.erase(det.TYuChannel.begin()+j);
						det.TYuNeighbour.erase(det.TYuNeighbour.begin()+j);
						det.TYuNo.erase(det.TYuNo.begin()+j);
						det.TYuRing.erase(det.TYuRing.begin()+j);
						det.TYuTheta.erase(det.TYuTheta.begin()+j);
						break;
					}
				}
			}
			if(det.TYuMul>=i) break;
		}
  */

		// Zero-degree DSSD, vertical strips
		for (Int_t i=0;i<NZdxChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NZdxChannels;j++){
        		if(Zdx[j] > maxE){
            		maxE = Zdx[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TZdxMul++;
    			det.TZdxEnergy.push_back(Zdx[maxCh]);
				det.TZdxChannel.push_back(maxCh);
				xpos = (maxCh-8)*3.+1.5;
				det.TZdxPos.push_back(xpos);
				//rndm = 0.99*fRandom.Rndm(); //random number between 0 and 0.99 for each event
				//theta = TMath::RadToDeg()*atan((geoM.SdInnerRadius*(24.-maxCh-rndm)+geoM.SdOuterRadius*(maxCh+rndm))/24./(geoM.Sd1Distance+14.8));
				//det.TSd2Theta.push_back(theta); //AS theta angle for Sd (24 - number of rings)

				Zdx[maxCh] = 0.;
			}
			else break;
    	}
						
		// Zero-degree DSSD, horizontal strips
		for (Int_t i=0;i<NZdyChannels;i++){
    		maxE = 0.;
    		maxCh = -1;
    		for (Int_t j=0;j<NZdyChannels;j++){
        		if(Zdy[j] > maxE){
            		maxE = Zdy[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TZdyMul++;
    			det.TZdyEnergy.push_back(Zdy[maxCh]);
				det.TZdyChannel.push_back(maxCh);
				ypos = (maxCh-8)*3.+1.5;
				det.TZdyPos.push_back(ypos);
				// phi = 180.-360.*maxCh/32.;
				// rndm = 0.99*fRandom.Rndm(); //random number between 0 and 0.99 for each event
				// det.TSd2Phi.push_back(phi-11.25*rndm);

				Zdy[maxCh] = 0.;
			}
			else break;
    	}
		
		det.TZdMul = (det.TZdyMul<det.TZdxMul) ? det.TZdyMul: det.TZdxMul;
		for(Int_t i=0; i<det.TZdMul; i++){
			radius = TMath::Sqrt(det.TZdxPos.at(i)*det.TZdxPos.at(i)+det.TZdyPos.at(i)*det.TZdyPos.at(i));
			det.TZdR.push_back(radius);
			theta = TMath::ATan2(radius,(geoM.Sd1Distance+29.3))*TMath::RadToDeg();
			det.TZdTheta.push_back(theta);
			phi = TMath::ATan2(det.TZdyPos.at(i),det.TZdxPos.at(i))*TMath::RadToDeg();
			det.TZdPhi.push_back(phi);
		}
		// CsI
		for (Int_t i=0;i<NCsIChannels;i++){
    		maxE = 0.;
    		maxCh = 0;
    		for (Int_t j=0;j<NCsIChannels;j++){
        		if(CsI1[j] > maxE){
            		maxE = CsI1[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TCsI1Mul++;
    			if (det.TYdMul>0){
	      			int m = (det.TYdChannel.at(0)%16)/(16/NCsI1Group);
	      			det.TCsI1Energy.push_back((CsI1[maxCh]-CsI1Ped[maxCh])*CsI1Gain[m][maxCh]); 
	      			det.TCsI1Channel.push_back(maxCh);
          //if(maxCh==5 && m==3){ //Sector 5, ring 3.
          //  printf("%lf  %lf  %lf  %lf\n",CsI1[maxCh],CsI1Ped[maxCh],CsI1Gain[m][maxCh],det.TCsI1Energy.back());
          //}
					phi = 90.+1.75-360.*maxCh/16.;
					rndm = 22.4*fRandom.Rndm();
					phi = phi-11.2+rndm;
					phi = (phi<-180.)? phi+360. : phi;
					det.TCsI1Phi.push_back(phi);
	    		}
				CsI1[maxCh] = 0.;
			}
			else break;
    	}

		for (Int_t i=0;i<NCsIChannels;i++){
    		maxE = 0.;
    		maxCh = 0;
    		for (Int_t j=0;j<NCsIChannels;j++){
        		if(CsI2[j] > maxE){
            		maxE = CsI2[j];
            		maxCh = j;
        		}
    		}		
    
   			if(maxE>0.){ 
				det.TCsI2Mul++;
    			if (CsI2[maxCh] < 3840. && det.TYdMul>0){
					int m = (det.TYdChannel.at(0)%16)/(16/NCsI2Group);
	        		det.TCsI2Energy.push_back((CsI2[maxCh]-CsI2Ped[maxCh])*CsI2Gain[m][maxCh]);
	        		det.TCsI2Channel.push_back(maxCh);
					phi = 90.+1.75-360.*maxCh/16.;
					rndm = 22.4*fRandom.Rndm();
					phi = phi-11.2+rndm;
					phi = (phi<-180.)? phi+360. : phi;
					det.TCsI2Phi.push_back(phi);
	      		}
				CsI2[maxCh] = 0.;
			}
			else break;
    	}

  /*
		// Re-sorting Yd if hits don't match with CsI
		if(det.TCsI1Channel.size()>0&&det.TYdNo.size()>1&&calMesy.boolCsI1==true) // only check if YY1 has more than one hit and CsI has a hit and has been calibrated
		{
 			if (int(det.TCsI1Channel.at(0)/2)-det.TYdNo.at(0)!=0 && int(det.TCsI1Channel.at(0)/2)-det.TYdNo.at(1)==0)//checking if the CsI hit is behind the first or second hit in  YY1
			{
				std::swap(det.TYdEnergy.at(0),det.TYdEnergy.at(1));
				std::swap(det.TYdChannel.at(0),det.TYdChannel.at(1));
				std::swap(det.TYdNo.at(0),det.TYdNo.at(1));
				std::swap(det.TYdRing.at(0),det.TYdRing.at(1));
				std::swap(det.TYdTheta.at(0),det.TYdTheta.at(1));
			}
		}
  */

		// IC
		/*
		maxE=0; 
		maxCh = -1;
    	for (int i =0; i< NICChannels;i++) {
    		if (maxE<IC[i]){
      			maxE=IC[i];
      			maxCh = i;
			}
    	} //for
    	if(maxE>0.)
	*/
		for(int i=0; i< NICChannels; i++)
		{
			if(IC[i]>0.){
			//det.TICEnergy.push_back(maxE); //for filling the tree
			//det.TICChannel.push_back(maxCh);
			det.TICEnergy.push_back(IC[i]); //for filling the tree
			det.TICChannel.push_back(i);
			}
		}
		
		// SSB
		det.TSSBEnergy = SSBEnergy;

		// Scintillator
		det.TScEnergy = ScEnergy;
		
		if(TrEnergy[0]>0){
			for(int i=0; i<NTrChannels; i++){
				det.TTrEnergy.push_back(TrEnergy[i]);
			}
		}
		
 		*pdet = det;
	} //last bank
}

//---------------------------------------------------------------------------------
void HandleBOR_Mesytec(int run, int gFileNumber, int time, IDet* pdet, std::string CalibFile)
{
	if(CalibFile=="") printf("No calibration file specified!\n\n");
	calMesy.Load(CalibFile);
	calMesy.Print();

	geoM.ReadGeometry(calMesy.fileGeometry.data());
// ************************************************************************************

	treeFile->cd();
	// create a TTree
	if(gFileNumber==0){
	   	tree = new TTree("Iris","iris data");
		tree->Branch("det","IDet",pdet,32000,99);
	}
	else{
	   	tree = (TTree*)treeFile->Get("Iris");
		tree->SetBranchAddress("det",&pdet);
	}


// Temporary variables for calibration 
 	Int_t Chan=-1;
	double a,b;
	int g; //for ringwise calibration of CsI

//************** Calibrate IC, not yet implemented! *********************************
	FILE * pFile;
	FILE * logFile;
 	char buffer[32];

	// logfile
   	logFile = fopen("treeIris.log","w");
		
//*************** Calibrate IC ********************************
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

		for (int i =0;i<NICChannels;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
			ICPed[Chan] = a;
			ICGain[Chan] = b;
			printf("ICPed %lf ICgain %lf\n",ICPed[Chan],ICGain[Chan]);
     	}
     	fclose (pFile);
		printf("\n");
 	}
   	

//*************** Calibrate TRIFIC ****************************
	pFile = fopen(calMesy.fileTr.data(), "r");

	if (pFile == NULL || calMesy.boolTr==false) {
		//perror ("No file");
		fprintf(logFile,"No calibration file for TRIFIC. Skipping TRIFIC calibration.\n");
		printf("No calibration file for TRIFIC. Skipping TRIFIC calibration.\n");
		for (int i =0;i<NTrChannels;i++  ){
			TrPed[i] = 0.;
			TrGain[i] = 1.;
     	}
   	}  
 	else  {
		printf("Reading TRIFIC config file '%s'\n",calMesy.fileTr.data());
		// Skip first line
  		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		for (int i =0;i<16;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
			TrPed[Chan] = a;
			TrGain[Chan] = b;
			printf("TrPed %lf Trgain %lf\n",TrPed[Chan],TrGain[Chan]);
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
				printf("Sd2rOffset %lf Sd2rgain %lf\n",Sd2rOffset[Chan-64],Sd2rGain[Chan-64]);
			}
       		else if (usePeds){
				Sd2rPed[Chan-64] = a;
				Sd2rGain[Chan-64] =  b;  
				printf("Sd2rPed %lf Sd2rgain %lf\n",Sd2rPed[Chan-64],Sd2rGain[Chan-64]);
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
				printf("Sd2sOffset %lf Sd2sgain %lf\n",Sd2sOffset[Chan-96],Sd2sGain[Chan-96]);
			}
       		else if (usePeds){
				Sd2sPed[Chan-96] = a;
				Sd2sGain[Chan-96] = b;   
				printf("Sd2sPed %lf Sd2sgain %lf\n",Sd2sPed[Chan-96],Sd2sGain[Chan-96]);
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
				printf("Sd1sOffset %lf Sd1sgain %lf\n",Sd1sOffset[Chan-160],Sd1sGain[Chan-160]);
			}
       		else if (usePeds){
				Sd1sPed[Chan-160] = a;
				Sd1sGain[Chan-160] = b; 
				printf("Sd1sPed %lf Sd1sgain %lf\n",Sd1sPed[Chan-160],Sd1sGain[Chan-160]);
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
			SurOffset[i] = 0.;
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
       		if (!usePeds){
				SurOffset[Chan-320] = a;
				SurGain[Chan-320] = b; 
				printf("SurOffset %lf Surgain %lf\n",SurOffset[Chan-320],SurGain[Chan-320]);
			}
       		else if (usePeds){
				SurPed[Chan-320] = a;
				SurGain[Chan-320] = b;
				printf("SurPed %lf Surgain %lf\n",SurPed[Chan-320],SurGain[Chan-320]);
			}
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
   		for (int i =0;i<32;i++  ){
			SusPed[i] = 0.;
			SusOffset[i] = 0.;
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
				SusGain[Chan-352] = b; 
				printf("SusOffset %lf Susgain %lf\n",SusOffset[Chan-352],SusGain[Chan-352]);
			}
       		else if (usePeds){
				SusPed[Chan-352] = a;
				SusGain[Chan-352] = b; 
				printf("SusPed %lf Susgain %lf\n",SusPed[Chan-352],SusGain[Chan-352]);
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
			YdPedestal[i] = 0.;
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
 			if(!usePeds){
			   	YdOffset[Chan-192] = a;
 				printf("gain %lf offset %lf\t",b,a);
			}
			else{
			   	YdPedestal[Chan-192] = a;
 				printf("gain %lf pedestal %lf\t",b,a);
			}
			if((i+1)%4==0) printf("\n");
    	}
    	fclose (pFile);
		printf("\n");
  	}
//************************************************************************
//**************** Calibrate Yu ****************************************
	Chan=-1;
   	pFile = fopen (calMesy.fileYu.data(), "r");

	if (pFile == NULL || calMesy.boolYu==false) {
		fprintf(logFile,"No calibration file for Yu. Skipping Yu calibration.\n");
		printf("No calibration file for Yu. Skipping Yu calibration.\n");
		for (int i =0;i<NYuChannels;i++  ){
			YuOffset[i] = 0.;
			YuPedestal[i] = 0.;
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
 			if(!usePeds){
			   	YuOffset[Chan-384] = a;
 				printf("gain %lf offset %lf\t",b,a);
			}
			else{
			   	YuPedestal[Chan-384] = a;
 				printf("gain %lf pedestal %lf\t",b,a);
			}
			if((i+1)%4==0) printf("\n");
    	}
		fclose (pFile);
		printf("\n");
  	}
//************************************************************************

//**************** Calibrate Zd horizontal strips ****************************************
	Chan=-1;

	pFile = fopen (calMesy.fileZdy.data(), "r");

   	if (pFile == NULL || calMesy.boolZdy==false) {
		fprintf(logFile,"No calibration file for Zd horizontal strips. Skipping Zdy calibration.\n");
		printf("No calibration file for Zd horizontal strips. Skipping Zdy calibration.\n");
   		for (int i =0;i<16;i++  ){
			ZdyPed[i] = 0.;
			ZdyOffset[i] = 0.;
			ZdyGain[i] = 1.;  
		}
	}  
 	else {
		printf("Reading Zdy config file '%s'\n",calMesy.fileZdy.data());
		// Skip first line
  		fscanf(pFile,"%s",buffer);
  		fscanf(pFile,"%s",buffer);
 		fscanf(pFile,"%s",buffer);

		for (int i =0;i<16;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if(!usePeds){
				ZdyOffset[Chan-64] = a;
				ZdyGain[Chan-64] =  b;  
				printf("%d ZdyOffset %lf Zdygain %lf\n",Chan-64,ZdyOffset[Chan-64],ZdyGain[Chan-64]);
			}
       		else if (usePeds){
				ZdyPed[Chan-64] = a;
				ZdyGain[Chan-64] =  b;  
				printf("%d ZdyPed %lf Zdygain %lf\n",Chan-64,ZdyPed[Chan-64],ZdyGain[Chan-64]);
     		}
		}
     	fclose (pFile);
		printf("\n");
 	}
//************************************************************************

//**************** Calibrate Zd vertical strips ****************************************
	Chan=-1;
   	pFile = fopen (calMesy.fileZdx.data(), "r");

 	if (pFile == NULL || calMesy.boolZdx==false) {
		fprintf(logFile,"No calibration file for Zd vertical strips. Skipping Zdx calibration.\n");
		printf("No calibration file for Zd vertical strips. Skipping Zdx calibration.\n");
   		for (int i =0;i<16;i++  ){
			ZdxPed[i] = 0.;
			ZdxOffset[i] = 0.;
			ZdxGain[i] = 1.;  
		}
	}  
 	else  {
		printf("Reading Zdx config file '%s'\n",calMesy.fileZdx.data());
		// Skip first line
	  	fscanf(pFile,"%s",buffer);
	  	fscanf(pFile,"%s",buffer);
	 	fscanf(pFile,"%s",buffer);

		for (int i =0;i<16;i++  ){
       		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
       		if (!usePeds){
				ZdxOffset[Chan-96] = a;
				ZdxGain[Chan-96] = b;   
				printf("%d ZdxOffset %lf Zdxgain %lf\n",Chan-64,ZdxOffset[Chan-96],ZdxGain[Chan-96]);
			}
       		else if (usePeds){
				ZdxPed[Chan-96] = a;
				ZdxGain[Chan-96] = b;   
				printf("%d ZdxPed %lf Zdxgain %lf\n",Chan-96,ZdxPed[Chan-96],ZdxGain[Chan-96]);
			}
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

	printf("Finished HandleBOR_Mesytec\n");
//************************************************************************


}

void HandleEOR_Mesytec(int run, int time)
{
  printf(" in Mesytec EOR\n");
}
