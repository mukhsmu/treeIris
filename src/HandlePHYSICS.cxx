// HandlePhysics.cxx
// calculating Physics variables requiring data from more than one detector

#include <iostream>
#include <assert.h>
#include <fstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TMath.h>
#include "TCutG.h"
#include "HandleMesytec.h"
#include "HandlePHYSICS.h"
#include "Globals.h"
#include "eloss.h"
#include "nucleus.h"
#include "runDepPar.h"
#include "CalibPHYSICS.h"
#include "Graphsdedx.h"
#include "geometry.h"

//#define pd // (p,d) reaction analysis
extern TEvent *IrisEvent;
extern TFile* treeFile;
extern TTree* tree;

//#define csv_output //output a selected dataset into a csv file, for analysis with a different program, such as python
#ifdef csv_output
ofstream csv_file;
#endif

CalibPHYSICS calPhys;
geometry geoP;
Graphsdedx dedx_p, dedx_d, dedx_t, dedx_i[3];

float ICELossCorr;

const int Nchannels = 24;
const int binlimit = 1900;
 	
Double_t eATgt[3][100], eAIso[3][100], eAWndw[3][100], eAAg[3][100];	
Double_t dedxATgt[3][100], dedxAIso[3][100], dedxAWndw[3][100], dedxAAg[3][100];	
Double_t eBSi[3][100], eBTgt[3][100], eBSiO2[3][100], eBB[3][100], eBP[3][100], eBAl[3][100], eBMy[3][100], eBCsI[3][100];	
Double_t dedxBSi[3][100], dedxBTgt[3][100], dedxBSiO2[3][100], dedxBB[3][100], dedxBP[3][100], dedxBAl[3][100], dedxBMy[3][100], dedxBCsI[3][100];	
Double_t ebTgt[3][100], ebSi[3][100], ebAl[3][100], ebB[3][100], ebMy[3][100], ebP[3][100], ebCsI[3][100], ebSiO2[3][100];	
Double_t dedxbTgt[3][100], dedxbSi[3][100], dedxbAl[3][100], dedxbB[3][100], dedxbMy[3][100], dedxbP[3][100], dedxbCsI[3][100], dedxbSiO2[3][100];	

double EBAC = 0.; //Beam energy from accelerator
runDep runDepPar[3]; // run dependant parameters
nucleus beam[3]; // beam particle properties
nucleus target; // target properties
nucleus lej; // light ejetcile properties
nucleus hej[3]; // heavy ejectile properties  
Int_t useYCalc = 0;//Use calculated YY1 energy  : do we need to ? :Jaspreet
Double_t PResid; //Momentum of residue
Double_t PBeam; // Calculated beam momentum after scattering off Ag
Double_t PA; //Beam momentum before reaction
Double_t Pb1; //Light ejectile momentum
Double_t Pb2; //Light ejectile momentum

Double_t A,B,C; //Used for quadratic equations
Double_t MBeam = 0.; // Beam mass
Double_t kBF = 108.904752/(MBeam/931.494013); //Ratio of Beam particle mass and 109-Ag foil nucleus mass
const Double_t MFoil = 931.494013*108.; //Ag foil mass AS
//Double_t geoP.FoilThickness = 5.7;                      /// ?????????? confrim if 5.7 or 5.44 --Jaspreet
Double_t energy = 0;
Double_t cosTheta = 1.;// S3 angle cosine
Double_t EBeam;
Double_t betaCM, gammaCM; //CM velocity

Double_t mA= 0.; //11.;//Beam mass //Reassigned in HandleBOR
Double_t ma = 0.; //target mass
Double_t mb= 0.; //1.;//Light ejectile mass
Double_t mB= 0.; // Heavy ejectile mass

Double_t b1,j,thetaD,E1,Eb1,Eb2,thetaR,ECsI1,ECsI2,EYY1;
Double_t Q1,Q2,thetaCM1,thetaCM2;
Double_t Pb1y,Pb2y,Pb1xcm,Pb2xcm;
	
//TCutG * protons;
//TCutG * deuterons;

int nGate = 0; // for selecting incoming ion

//---------------------------------------------------------------------------------
void HandleBOR_PHYSICS(int run, int time, IDet *det, TString CalibFile)
{
  	printf("In HandleBOR_PHYSICS...\n");	
	calPhys.Load(CalibFile);
	calPhys.Print();

//	TFile *fgates = new TFile(calPhys.fileGates.data());
//   	printf("opened file %s\n",calPhys.fileGates.data());
//  
//	protons = (TCutG*)fgates->FindObjectAny("proton");
//  	if(!protons) printf("No proton gate.\n");  
// 	deuterons = (TCutG*)fgates->FindObjectAny("deuteron");   // to be define for Mg experiment using 20Ne runs
//  	if(!deuterons) printf("No deuteron gate.\n");  
//	else deuterons->SetName("deuterons");
//	printf("Grabbed gates.\n");
	
	#ifdef csv_output  
	  	csv_file.open("S1147pd.csv", ios::out);
	  	assert(csv_file.is_open());
	#endif

	geoP.ReadGeometry(calPhys.fileGeometry.data());
	geoP.Print();
	
	// Time dependent correction of IC energy loss 
	FILE * pFile;
	int Chan = 0;
	double a,b;
	int run_for_corr = 0;
	
	if (pFile == NULL || calPhys.boolTCorrIC==false) {
		//fprintf(logFile,"No time dependent correction for IC energy loss. Skipping correction.\n");
		printf("No time dependent correction for IC energy loss. Skipping correction.\n");
		ICELossCorr =1.;
	}
	else  {
		printf("Reading config file '%s'\n",calPhys.fileTCorrIC.data());

		while (!feof(pFile)){
    		fscanf(pFile,"%d%lf%lf",&Chan,&a,&b);
			if(Chan==run){ 
				run_for_corr = Chan;
				ICELossCorr = b; 
			}
    	}
    	fclose (pFile);	

		if(run_for_corr==0){ 
			printf("Run %d not in list. No correction applied!\n",run);
		}
		else{
			printf("Run: %d\tIC Gain correction: %f\n\n",Chan,ICELossCorr);
		}
  	}

	for(int i=0; i<3; i++){	
		if(calPhys.boolRunDepPar[i]){
			runDepPar[i].setRunDepPar(calPhys.fileRunDepPar[i]);// setting run dependent parameters.
			runDepPar[i].Print();

			beam[i].getInfo(runDepPar[i].nA);
			target.getInfo(runDepPar[i].na);
			hej[i].getInfo(runDepPar[i].nB);
			lej.getInfo(runDepPar[i].nb);

			EBAC = runDepPar[i].EBAC;
			mA = beam[i].mass; //Beam mass //Reassigned in HandleBOR
			ma = target.mass;
			mb = lej.mass; //Light ejectile mass
			mB = hej[i].mass;
      		printf("mA=%f, ma=%f, mb=%f, mB=%f\n",mA,ma,mb,mB); 
			kBF = MFoil/mA;
	
			printf("Beam energy: %f\n", runDepPar[i].energy);
			// printf("Target thickness: %f\n",geoP.TargetThickness);
			
			if(calPhys.boolIdedx[i]==kTRUE){
				printf("\n\nLoading dedx Graphs for ion %d...\n",(i+1));
				dedx_i[i].Load(calPhys.fileIdedx[i]);
				dedx_i[i].Print();
				if(dedx_i[i].boolAg==kTRUE) loadELoss(dedx_i[i].Ag,eAAg[i],dedxAAg[i],mA);	
				if(dedx_i[i].boolTgt==kTRUE) loadELoss(dedx_i[i].Tgt,eATgt[i],dedxATgt[i],mA);	
				if(dedx_i[i].boolSi==kTRUE) loadELoss(dedx_i[i].Si,eBSi[i],dedxBSi[i],mB);	
				if(dedx_i[i].boolAl==kTRUE) loadELoss(dedx_i[i].Al,eBAl[i],dedxBAl[i],mB);	
				if(dedx_i[i].boolB==kTRUE) loadELoss(dedx_i[i].B, eBB[i],dedxBB[i],mB);	
				if(dedx_i[i].boolP==kTRUE) loadELoss(dedx_i[i].P, eBP[i],dedxBP[i],mB);	
				if(dedx_i[i].boolSiO2==kTRUE) loadELoss(dedx_i[i].SiO2,eBSiO2[i],dedxBSiO2[i],mB);	
				if(dedx_i[i].boolIso==kTRUE) loadELoss(dedx_i[i].Iso,eAIso[i],dedxAIso[i],mB);	
				if(dedx_i[i].boolWndw==kTRUE) loadELoss(dedx_i[i].Wndw,eAWndw[i],dedxAWndw[i],mB);	
			}
			//Needs to be moved!
			const Double_t ICLength=22.9*0.062; //cm*mg/cm^3 at 19.5 Torr
			const Double_t ICWindow1=0.03*3.44*0.1; //mu*g/cm^3*0.1
			const Double_t ICWindow2=0.05*3.44*0.1; //mu*g/cm^3*0.1
	
			if (eAIso[i]&&eAWndw[i]){
				Double_t temp_energy = runDepPar[i].energy;
				runDepPar[i].energy = runDepPar[i].energy-eloss(runDepPar[i].energy,ICWindow1,eAWndw[i],dedxAWndw[i]);  
				runDepPar[i].energy = runDepPar[i].energy-eloss(runDepPar[i].energy,ICLength,eAIso[i],dedxAIso[i])/ICELossCorr;  
				runDepPar[i].energy = runDepPar[i].energy-eloss(runDepPar[i].energy,ICWindow2,eAWndw[i],dedxAWndw[i]);  
				printf("Energy loss in IC (including windows): %f\n" ,temp_energy-runDepPar[i].energy);
			}
			else printf("Energy loss in IC/IC windows not specified. EBeam not corrected.");


			if (eAAg[i]){
				Double_t temp_energy = runDepPar[i].energy;
				runDepPar[i].energy = runDepPar[i].energy-eloss(runDepPar[i].energy,geoP.FoilThickness,eAAg[i],dedxAAg[i]);  
				printf("Energy loss in silver foil: %f\n" ,temp_energy-runDepPar[i].energy);
			}
			else printf("Energy loss in silver foil not specified. EBeam not corrected.");

			printf("Energy after silver foil: %f\n",runDepPar[i].energy);

			if (eATgt[i]){
				Double_t temp_energy = runDepPar[i].energy;
				runDepPar[i].energy = runDepPar[i].energy-eloss(runDepPar[i].energy,geoP.TargetThickness/2.,eATgt[i],dedxATgt[i]);  
				printf("Energy loss in half target: %f\n" ,temp_energy-runDepPar[i].energy);
			}
			else printf("Energy loss in target not specified. EBeam not corrected.");

			runDepPar[i].momentum = sqrt(runDepPar[i].energy*runDepPar[i].energy+2.*runDepPar[i].energy*beam[i].mass);//beam momentum
			runDepPar[i].beta = runDepPar[i].momentum/(runDepPar[i].energy + beam[i].mass + target.mass);
			runDepPar[i].gamma = 1./sqrt(1.-runDepPar[i].beta*runDepPar[i].beta);

			EBeam= runDepPar[i].energy;
			PA = runDepPar[i].momentum;//beam momentum
			betaCM = runDepPar[i].beta;
			gammaCM = runDepPar[i].gamma;
	
			printf("Resulting energy: %f\n",EBeam);
			printf("Beam momentum: %f\n",PA);
			printf("Beam beta: %f\tBeam gamma: %f\n",betaCM,gammaCM);
			printf("MBeam: %f\t MFoil: %f\t kBF: %f\n",beam[i].mass,MFoil,kBF);
			printf("beam mass: %f MeV (%f)\ttarget mass: %f MeV (%f)\n",mA,mA/931.494061,ma,ma/931.494061);
			printf("heavy ejectile mass: %f MeV (%f)\tlight ejectile mass: %f MeV (%f)\n",mB,mB/931.494061,mb,mb/931.494061);
		}
	}
	//---------------for protons in various layers---Jaspreet ----------------------------------
	if(calPhys.boolPdedx==kTRUE){
		printf("\n\nLoading dedx Graphs for proton...\n");
		dedx_p.Load(calPhys.filePdedx);
		dedx_p.Print();
		if(dedx_p.boolSi==kTRUE) loadELoss(dedx_p.Si,ebSi[0],dedxbSi[0],mb);	
		if(dedx_p.boolAl==kTRUE) loadELoss(dedx_p.Al,ebAl[0],dedxbAl[0],mb);	
		if(dedx_p.boolB==kTRUE) loadELoss(dedx_p.B,ebB[0],dedxbB[0],mb);	
		if(dedx_p.boolP==kTRUE) loadELoss(dedx_p.P,ebP[0],dedxbP[0],mb);	
		if(dedx_p.boolMy==kTRUE) loadELoss(dedx_p.My,ebMy[0],dedxbMy[0],mb);	
		if(dedx_p.boolTgt==kTRUE) loadELoss(dedx_p.Tgt,ebTgt[0],dedxbTgt[0],mb);	
	}
	//---------------for deuterons in various layers---Jaspreet ----------------------------------
	if(calPhys.boolDdedx==kTRUE){
		printf("\n\nLoading dedx Graphs for deuteron...\n");
		dedx_d.Load(calPhys.fileDdedx);
		dedx_d.Print();
		if(dedx_d.boolSi==kTRUE) loadELoss(dedx_d.Si,ebSi[1],dedxbSi[1],mb);	
		if(dedx_d.boolAl==kTRUE) loadELoss(dedx_d.Al,ebAl[1],dedxbAl[1],mb);	
		if(dedx_d.boolB==kTRUE)  loadELoss(dedx_d.B,ebB[1],dedxbB[1],mb);	
		if(dedx_d.boolP==kTRUE)  loadELoss(dedx_d.P,ebP[1],dedxbP[1],mb);	
		if(dedx_d.boolMy==kTRUE) loadELoss(dedx_d.My,ebMy[1],dedxbMy[1],mb);	
		if(dedx_d.boolTgt==kTRUE) loadELoss(dedx_d.Tgt,ebTgt[1],dedxbTgt[1],mb);	
	}

//---------------for tritons in various layers---Jaspreet ----------------------------------
	if(calPhys.boolTdedx==kTRUE){
		printf("\n\nLoading dedx Graphs for deuteron...\n");
		dedx_t.Load(calPhys.fileDdedx);
		dedx_t.Print();
		if(dedx_t.boolSi==kTRUE) loadELoss(dedx_t.Si,ebSi[2],dedxbSi[2],mb);	
		if(dedx_t.boolAl==kTRUE) loadELoss(dedx_t.Al,ebAl[2],dedxbAl[2],mb);	
		if(dedx_t.boolB==kTRUE)  loadELoss(dedx_t.B,ebB[2],dedxbB[2],mb);	
		if(dedx_t.boolP==kTRUE)  loadELoss(dedx_t.P,ebP[2],dedxbP[2],mb);	
		if(dedx_t.boolMy==kTRUE) loadELoss(dedx_t.My,ebMy[2],dedxbMy[2],mb);	
		if(dedx_t.boolTgt==kTRUE) loadELoss(dedx_t.Tgt,ebTgt[2],dedxbTgt[2],mb);	
	}
//--------------------------------------------------------------------------------
		printf("End of HandleBOR_Physics\n");
} //HandleBOR_Physics

void HandlePHYSICS(IDet *det)
{
	if(calPhys.boolEssential==false || det->TSd1rMul==0){ // not enough calibration information or no hit in first S3
	   tree->Fill();
   	   return;
	} 
  	
	cosTheta = cos(TMath::DegToRad()* (det->TSdTheta.at(0)));

	// Selecting an incoming isotope:	
 	if (calPhys.boolICGates==kFALSE) nGate=0;   
 	else if (calPhys.boolICGates==kTRUE) {  
		for(int i=0; i<3; i++){  
 			if ((det->TICEnergy>runDepPar[i].ICmin) && (det->TICEnergy<runDepPar[i].ICmax)) { 
				nGate = i; 
  			}
		}
	} 

	EBAC = runDepPar[nGate].EBAC;
	mA = beam[nGate].mass; 
	ma = target.mass; 
	mB = hej[nGate].mass;
	mb = lej.mass;
	kBF = MFoil/mA;
	EBeam= runDepPar[nGate].energy;
	betaCM = runDepPar[nGate].beta;
	gammaCM = runDepPar[nGate].gamma;
	PA = runDepPar[nGate].momentum; // beam momentum
	
 	//adding dead layer energy losses
	//Sd2 ring side
	if(det->TSd1rEnergy.size()>0 && det->TSd2rEnergy.size()>0){
		if(det->TSd1rEnergy.at(0)>0. && det->TSd2rEnergy.at(0)>0.){
			energy = det->TSd2rEnergy.at(0);
			if (eBB[nGate]) energy = energy+elossFi(energy,0.1*2.35*0.5/cosTheta,eBB[nGate],dedxBB[nGate]); //boron junction implant
			if (eBAl[nGate]) energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,eBAl[nGate],dedxBAl[nGate]); //first metal
			if (eBSiO2[nGate]) energy = energy+elossFi(energy,0.1*2.65*2.5/cosTheta,eBSiO2[nGate],dedxBSiO2[nGate]); //SiO2
			if (eBAl[nGate]) energy = energy+elossFi(energy,0.1*2.7*1.5/cosTheta,eBAl[nGate],dedxBAl[nGate]); //second metal
			//Sd1 ring side
			if (eBAl[nGate]) energy = energy+elossFi(energy,0.1*2.7*1.5/cosTheta,eBAl[nGate],dedxBAl[nGate]); //second metal
			if (eBSiO2[nGate]) energy = energy+elossFi(energy,0.1*2.65*2.5/cosTheta,eBSiO2[nGate],dedxBSiO2[nGate]); //SiO2
			if (eBAl[nGate]) energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,eBAl[nGate],dedxBAl[nGate]); //first metal
			energy = energy + det->TSd1rEnergy.at(0);// energy lost and measured in Sd1

			//sector side
			if (eBP[nGate]) energy = energy+elossFi(energy,0.1*1.822*0.5/cosTheta,eBP[nGate],dedxBP[nGate]); //phosphorus implant
			if (eBAl[nGate]) det->TSdETot = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,eBAl[nGate],dedxBAl[nGate]); //metal
		}
	}

	PResid = sqrt(2.*det->TSdETot*mA);     //Beam momentum in MeV/c
	A = kBF-1.;                              //Quadratic equation parameters
    B = 2.0*PResid* cos(TMath::DegToRad()*det->TSdTheta.at(0));
    C = -1.*(kBF+1)*PResid*PResid; 
    if (A!=0)    PBeam = (sqrt(B*B-4.*A*C)-B)/(2*A);
  	//to calculate residue energy from beam
    
	IrisEvent->fEB = PBeam*PBeam/(2.*mA);
	//det->TBE = PBeam*PBeam/(2.*mA);
   
    // printf("thetaCM: %f\n",det->TSdThetaCM);
    if (eAAg[nGate]) IrisEvent->fEB=  IrisEvent->fEB + elossFi(det->TSdETot,geoP.FoilThickness/2.,eAAg[nGate],dedxAAg[nGate]); //energy loss from the end of H2 to the center of Ag.
    det->TSdThetaCM = TMath::RadToDeg()*atan(tan(TMath::DegToRad()*det->TSdTheta.at(0))/sqrt(gammaCM-gammaCM*betaCM*(mA+IrisEvent->fEB)/(PBeam*cos(TMath::DegToRad()*det->TSdTheta.at(0)))));// check if this is still correct for H2 target tk
	if(det->TYdMul==0||det->TCsI2Mul==0){
	   tree->Fill();
   	   return;
	}
 // Calculate Q-value from YY1 and CsI// 
	//if ((deuterons->IsInside(det->TCsI2Energy.at(0),det->TYdEnergy.at(0)*cos(det->TYdTheta.at(0)*TMath::Pi()/180.))) && ((det->TYdEnergy.at(0)>0.2)  && (det->TCsI2Energy.at(0)>0.6 ))) {    //check if in the proton/deuteron gate
	if (det->TYdEnergy.size()>0&&det->TYdRing.size()>0) {    //check if in the proton/deuteron gate
	  	thetaR =( atan((geoP.YdInnerRadius+((det->TYdRing.at(0)+1)*(geoP.YdOuterRadius-geoP.YdInnerRadius)/16))/geoP.YdDistance) + atan((geoP.YdInnerRadius+((det->TYdRing.at(0))*(geoP.YdOuterRadius-geoP.YdInnerRadius)/16))/geoP.YdDistance) )/2.;
	  	thetaD = thetaR*TMath::RadToDeg();
		EYY1 = det->TYdEnergy.at(0);
	}	
	if (det->TCsI1Energy.size()>0&&det->TYdEnergy.size()>0&&det->TYdRing.size()>0) {    //check if in the proton/deuteron gate
	
	    ECsI1= det->TCsI1Energy.at(0);
	
	  	if (mb == target.mass) //proton energy loss in dead layers between YY1 and CsI                                                                                       
	    {
	      	if (ebMy[0])  ECsI1= ECsI1+elossFi(ECsI1,0.1*1.4*6./cos(thetaR),ebMy[0],dedxbMy[0]); //Mylar                                                                                  
	      	else std::cout << "ebMy doesn't exist"<< std::endl;
	      	if (ebAl[0])  ECsI1= ECsI1+elossFi(ECsI1,0.1*2.702*0.3/cos(thetaR),ebAl[0],dedxbAl[0]); //0.3 u Al                                                                            
	      	else std::cout << "ebAl doesn't exist"<< std::endl;
	      	if (ebP[0])  ECsI1= ECsI1+elossFi(ECsI1,0.1*1.88219*0.1/cos(thetaR),ebP[0],dedxbP[0]); // 0.1Phosphorus                                                                      
	      	else std::cout << "ebP doesn't exist"<< std::endl;
	    }
	
		Eb1= ECsI1+EYY1; //use measured Yd // change june28
			
	   	if (mb == target.mass){
	      	if (ebSi[0])  Eb1= Eb1+elossFi(Eb1,0.1*2.32*0.35/cos(thetaR),ebSi[0],dedxbSi[0]); //0.3 u Al + 1 um B equivalent in 0.35 um Si                                                            
	    	else std::cout << "ebSi doesn't exist"<< std::endl;
	    	if (ebTgt[0])  Eb1= Eb1+elossFi(Eb1,geoP.TargetThickness/2./cos(thetaR),ebTgt[0],dedxbTgt[0]); //deuteron energy  in mid target midtarget
	     	else std::cout << "ebD2 doesn't exist"<< std::endl;
		}
		Pb1 = sqrt(Eb1*Eb1+2.*Eb1*mb);
		Pb1y = Pb1*sin(thetaR);
		Pb1xcm = gammaCM*betaCM*(Eb1+mb)- gammaCM*Pb1*cos(thetaR);
	 	Q1 = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb1)+2.*PA*Pb1*cos(thetaR)+2.*(EBeam+mA+ma-Eb1-mb)*ma);  //Alisher's equation 
	  	IrisEvent->fQv1 = Q1;
		thetaCM1 = TMath::RadToDeg()*atan(Pb1y/Pb1xcm);
		thetaCM1 = (thetaCM1<0) ? thetaCM1+180. : thetaCM1;
		IrisEvent->fThetacm1 = thetaCM1;
	}
	
	if (det->TCsI2Energy.size()>0&&det->TYdEnergy.size()>0&&det->TYdRing.size()>0) {    //check if in the proton/deuteron gate
	
	    ECsI2= det->TCsI2Energy.at(0);
	
	  	if (mb == target.mass) //proton energy loss in dead layers between YY1 and CsI                                                                                       
	    {
	      	if (ebMy[0])  ECsI2= ECsI2+elossFi(ECsI2,0.1*1.4*6./cos(thetaR),ebMy[0],dedxbMy[0]); //Mylar                                                                                  
	      	else std::cout << "ebMy doesn't exist"<< std::endl;
	      	if (ebAl[0])  ECsI2= ECsI2+elossFi(ECsI2,0.1*2.702*0.3/cos(thetaR),ebAl[0],dedxbAl[0]); //0.3 u Al                                                                            
	      	else std::cout << "ebAl doesn't exist"<< std::endl;
	      	if (ebP[0])  ECsI2= ECsI2+elossFi(ECsI2,0.1*1.88219*0.1/cos(thetaR),ebP[0],dedxbP[0]); // 0.1Phosphorus                                                                      
	      	else std::cout << "ebP doesn't exist"<< std::endl;
	    }
	
		Eb2= ECsI2+EYY1; //use measured Yd // change june28
	
	   	if (mb == target.mass){
	      	if (ebSi[0])  Eb2= Eb2+elossFi(Eb2,0.1*2.32*0.35/cos(thetaR),ebSi[0],dedxbSi[0]); //0.3 u Al + 1 um B equivalent in 0.35 um Si                                                            
	    	else std::cout << "ebSi doesn't exist"<< std::endl;
	    	if (ebTgt[0])  Eb2= Eb2+elossFi(Eb2,geoP.TargetThickness/2./cos(thetaR),ebTgt[0],dedxbTgt[0]); //deuteron energy  in mid target midtarget
	     	else std::cout << "ebD2 doesn't exist"<< std::endl;
		}
		Pb2 = sqrt(Eb2*Eb2+2.*Eb2*mb);
		Pb2y = Pb2*sin(thetaR);
		Pb2xcm = gammaCM*betaCM*(Eb2+mb)- gammaCM*Pb2*cos(thetaR);
	 	Q2 = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb2)+2.*PA*Pb2*cos(thetaR)+2.*(EBeam+mA+ma-Eb2-mb)*ma);  //Alisher's equation 
	  	IrisEvent->fQv2 = Q2;
		thetaCM2 = TMath::RadToDeg()*atan(Pb2y/Pb2xcm);
		thetaCM2 = (thetaCM2<0) ? thetaCM2+180. : thetaCM2;
		IrisEvent->fThetacm2 = thetaCM2;
	}
	
	//	My eqs  
	
	 	// cout<<"Beam energy used for QValue is "<< EBeam <<endl;
		//	Q = mA+ma-mb-sqrt((EBeam+mA+ma-Eb-mb)*(EBeam+mA+ma-Eb-mb)-(PA*PA+Pb*Pb-2.*PA*Pb*cos(thetaR))); //Equivalent to the previous equation
	
		// EB = EBeam+mA+ma-Eb-mb; // rel. kinetic energy of recoil.
		// PB = sqrt(PA*PA+Pb*Pb-2.*PA*Pb*cos(thetaR));// momentum of the recoil
		// mB = sqrt(EB^2-PB^2);
	  	//det->QValue1 = Q1;
	//	IrisEvent->fThetacm = TMath::RadToDeg()*atan(Pby/Pbxcm);
	//	if (IrisEvent->fThetacm < 0) 
	//  		IrisEvent->fThetacm = IrisEvent->fThetacm+180.;
	 		#ifdef csv_output
				csv_file << det->QValue <<"," << det->TYdEnergy<<","<<IrisEvent->fPart.fThetacm <<","<< det->TYdTheta<< endl;
			#endif
//	}//isinside proton (deuteron) gate

	tree->Fill();
}

void HandleEOR_PHYSICS(int run, int time)
{
 #ifdef csv_output
  csv_file.close();
#endif
  printf(" in Physics EOR\n");
}
 
