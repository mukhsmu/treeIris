//
// ROOT analyzer
//
// PHYSICS detector handling
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

#include <iostream>
#include <assert.h>
#include <fstream>

#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include "HandleMesytec.h"
#include "HandlePHYSICS.h"
#include "Globals.h"
#include "eloss.h"
#include "runDepPar.h"
#include "TCutG.h"
//#include "TGraph.h"
#include "CalibPHYSICS.h"
#include "Graphsdedx.h"

//#define pd // (p,d) reaction analysis
extern TEvent *IrisEvent;

extern TFile* treeFile;

extern TTree* tree;

//#define csv_output //output a selected dataset into a csv file, for analysis with a different program, such as python

#ifdef csv_output
ofstream csv_file;
#endif

CalibPHYSICS calPhys;
Graphsdedx dedx_p, dedx_d, dedx_i[3];

const int Nchannels = 24;
const int binlimit = 1900;
//TH1D * hPhysics_Q[Nchannels] = {NULL}; // Q

TGraph *grIAg[3] ={NULL}; // stopping power in Ag foil
TGraph *grIH2[3] ={NULL}; // stopping in H2
TGraph *grISi[3] = {NULL};
TGraph *grIAl[3] = {NULL};
TGraph *grIB[3] = {NULL};
TGraph *grIP[3] = {NULL};
TGraph *grISiO2[3] = {NULL};

TGraph *grdSi = {NULL};
TGraph *grdAl = {NULL};
TGraph *grdB = {NULL};
TGraph *grdP = {NULL};
TGraph *grdMy = {NULL};
TGraph *grdD2 = {NULL};
TGraph *grdH2 = {NULL};

TGraph *grpSi = {NULL};
TGraph *grpAl = {NULL};
TGraph *grpB = {NULL};
TGraph *grpP = {NULL};
TGraph *grpMy = {NULL};
TGraph *grpD2 = {NULL};
TGraph *grpH2 = {NULL};

double EBAC = 0.; //Beam energy from accelerator
beam_t beam[3]; // beam properties such as, mass, charge and energy
target_t tar; //target properties
ejectile_t lej, hej; //ejectile properties
Int_t useYCalc = 1;//Use calculated YY1 energy  : do we need to ? :Jaspreet
//Double_t YdThickness[8] = {112.,109.*1.04/1.2,110.*1.5/1.59,106.*1.04/1.14,101.,109.*1.17/1.28,111.*1.12/1.22,103.};//Thicknesses for YY1 detectors 
Double_t YdThickness[8]= {104.65, 101.15, 106.125, 101.75,100.052, 105.65,102.48, 105.84}; // Nov 25,2014 // should probably go to parameter file (geometry.txt) 
Double_t PResid; //Momentum of residue
Double_t PBeam; // Calculated beam momentum after scattering off Ag
Double_t PA; //Beam momentum before reaction
Double_t Pb; //Light ejectile momentum

Double_t A,B,C; //Used for quadratic equations
Double_t MBeam = 0.; // Beam mass
Double_t kBF = 108.904752/(MBeam/kAmu); //Ratio of Beam particle mass and 109-Ag foil nucleus mass
const Double_t MFoil = 931.494013*108.; //Ag foil mass AS
Double_t foilTh = 5.7;                      /// ?????????? confrim if 5.7 or 5.44 --Jaspreet
Double_t energy = 0;
Double_t cosTheta = 1.;// S3 angle cosine
Double_t EBeam;
Double_t betaCM, gammaCM; //CM velocity
Double_t Sd1Thickness = 61.*2.321/10.; //Sd1 thickness ==> param file

Double_t mA= 0.; //11.;//Beam mass //Reassigned in HandleBOR
Double_t ma = 0.; //target mass
Double_t mb= 0.; //1.;//Light ejectile mass
Double_t mB= 0.; // Heavy ejectile mass
//Double_t a1, a2;
Double_t  b1,j,thetaD, Q,E1,Eb,thetaR,ECsI,ECsI1,EYY1;

TCutG * protons;
TCutG * deuterons;
//TGraph *Li11pp;//elastic scattering kinematics graph

float YdDistanceP = 0.; // distance from target in mm
float Yd1rP= 0., Yd2rP = 0. ; // inner and outer radii in mm
float Sd1DistanceP = 0., Sd2DistanceP = 0.; //distance from target in mm
float Sdr1P = 0., Sdr2P= 0.; //AS Inner and outer radii of an S3 detector (in mm).

int nGate = 0; // for selecting incoming ion

//---------------------------------------------------------------------------------
void HandleBOR_PHYSICS(int run, int time, IDet *det, TString CalibFile)
{
  	printf("In HandleBOR_PHYSICS...\n");	
	if(CalibFile=="0") printf("No calibration file specified");
	calPhys.Load(CalibFile);
	calPhys.Print();

	TFile *fgates = new TFile(calPhys.fileGates);
   	printf("opened file %s\n",calPhys.fileGates.Data());
  
	protons = (TCutG*)fgates->FindObjectAny("proton");
  	if(!protons) printf("No proton gate.\n");  
 	deuterons = (TCutG*)fgates->FindObjectAny("deuteron");   // to be define for Mg experiment using 20Ne runs
  	if(!deuterons) printf("No deuteron gate.\n");  
  	deuterons->SetName("deuterons");
	printf("Grabbed gates.\n");
	
	printf("Calculated dedx graphs in %s\n",calPhys.fileELoss.Data());
   	TFile *f1 = new TFile(calPhys.fileELoss);
 
#ifdef csv_output  
  	csv_file.open("S1147pd.csv", ios::out);
  	assert(csv_file.is_open());
#endif

//----------------stable beam for Mg experiment  i.e. 20Ne--------------------------------
	for(int i=0; i<3; i++){
		if(calPhys.boolIdedx[i]==kTRUE){
			printf("\n\nLoading dedx Graphs for ion %d...\n",(i+1));
			dedx_i[i].Load(calPhys.fileIdedx[i]);
			dedx_i[i].Print();
    		grIAg[i] = (TGraph*)f1->FindObjectAny(dedx_i[i].Ag);
			if(!grIAg[i]) printf("grIAg not found!\n");
    		grIH2[i] = (TGraph*)f1->FindObjectAny(dedx_i[i].H2);
			if(!grIH2[i]) printf("grIH2 not found!\n");
    		grISi[i] = (TGraph*)f1->FindObjectAny(dedx_i[i].Si);
			if(!grISi[i]) printf("grISi not found!\n");
    		grIAl[i] = (TGraph*)f1->FindObjectAny(dedx_i[i].Al); //density 2.7 g/cm^3 // i.5 um * 2.7 g/cm^3 = i.135 mg/cm^2 //i.8 um = i.216 mg/cm^2
			if(!grIAl[i]) printf("grIAl not found!\n");
    		grIB[i] = (TGraph*)f1->FindObjectAny(dedx_i[i].B);//i.5 um * 2.35 gm/cm^2 = i.1175 mg/cm^2                                              
			if(!grIB[i]) printf("grIB not found!\n");
    		grIP[i] = (TGraph*)f1->FindObjectAny(dedx_i[i].P); //i.5 um * 1.822 gm/cm^2 = i.i911 mg/cm^2
			if(!grIP[i]) printf("grIP not found!\n");
    		grISiO2[i] = (TGraph*)f1->FindObjectAny(dedx_i[i].SiO2);//2.5 um * 2.65 gm/cm^2 = i.6625 mg/cm^2
			if(!grISiO2[i]) printf("grISiO2 not found!\n");
		}
	}    
	
//---------------for protons in various layers---Jaspreet ----------------------------------
	if(calPhys.boolPdedx==kTRUE){
		printf("\n\nLoading dedx Graphs for proton...\n");
		dedx_p.Load(calPhys.filePdedx);
		dedx_p.Print();
 		grpSi = (TGraph*)f1->FindObjectAny(dedx_p.Si);
		if(!grpSi) printf("grpSi not found!\n");
   		grpAl = (TGraph*)f1->FindObjectAny(dedx_p.Al); //density 2.7 g/cm^3 // 0.5 um * 2.7 g/cm^3 = 0.135 mg/cm^2 //0.8 um = 0.216 mg/cm^2    
		if(!grpAl) printf("grpAl not found!\n");
    	grpB = (TGraph*)f1->FindObjectAny(dedx_p.B);//0.5 um * 2.35 gm/cm^2 = 0.1175 mg/cm^2                                                    
		if(!grpB) printf("grpB not found!\n");
    	grpP = (TGraph*)f1->FindObjectAny(dedx_p.P); //0.5 um * 1.822 gm/cm^2 = 0.0911 mg/cm^2                                                  
		if(!grpP) printf("grpP not found!\n");
  		grpMy = (TGraph*)f1->FindObjectAny(dedx_p.My); // 6 um * 1.4 gm/cm^2                                                                        
		if(!grpMy) printf("grpMy not found!\n");
		grpH2 = (TGraph*)f1->FindObjectAny(dedx_p.H2);
		if(!grpH2) printf("grpH2 not found!\n");
	}
//---------------for deuterons in various layers---Jaspreet ----------------------------------
	if(calPhys.boolDdedx==kTRUE){
		printf("\n\nLoading dedx Graphs for deuteron...\n");
		dedx_d.Load(calPhys.fileDdedx);
		dedx_d.Print();
		grdH2 = (TGraph*)f1->FindObjectAny(dedx_d.H2);
		if(!grdH2) printf("grdH2 not found!\n");
		grdAl = (TGraph*)f1->FindObjectAny(dedx_d.Al); //density 2.7 g/cm^3 // 0.5 um * 2.7 g/cm^3 = 0.135 mg/cm^2 //0.8 um = 0.216 mg/cm^2    
		if(!grdAl) printf("grdAl not found!\n");
	    grdB = (TGraph*)f1->FindObjectAny(dedx_d.B);//0.5 um * 2.35 gm/cm^2 = 0.1175 mg/cm^2                                                    
		if(!grdB) printf("grdB not found!\n");
	    grdP = (TGraph*)f1->FindObjectAny(dedx_d.P); //0.5 um * 1.822 gm/cm^2 = 0.0911 mg/cm^2                                                  
		if(!grdP) printf("grdP not found!\n");
	  	grdMy = (TGraph*)f1->FindObjectAny(dedx_d.My); // 6 um * 1.4 gm/cm^2                                                                        
		if(!grdMy) printf("grdMy not found!\n");
		grdD2 = (TGraph*)f1->FindObjectAny(dedx_d.D2);
		if(!grdD2) printf("grdD2 not found!\n");
		grdSi = (TGraph*)f1->FindObjectAny(dedx_d.Si);
		if(!grdSi) printf("grdSi not found!\n");
	}
//--------------------------------------------------------------------------------
	 
	f1->Close();
	printf("Closed f1\n");
	
	FILE * pFile;
	FILE * pwFile;
 	char buffer[32];
	
	pFile=fopen(calPhys.fileGeometry,"r");

	if (pFile == NULL) {
		perror ("Error opening file");
		fprintf(pwFile,"Error opening file");
   	}	
	printf("Reading config file '%s'\n",calPhys.fileGeometry.Data());
	
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
		if (strcmp(buffer,"YDD")==0)	YdDistanceP = v;
		if (strcmp(buffer,"YDR1")==0)	Yd1rP = v;
		if (strcmp(buffer,"YDR2")==0)	Yd2rP = v;
		if (strcmp(buffer,"SD1D")==0)	Sd1DistanceP = v;
		if (strcmp(buffer,"SD2D")==0)	Sd2DistanceP = v;
		if (strcmp(buffer,"SDR1")==0)	Sdr1P = v;
		if (strcmp(buffer,"SDR2")==0)	Sdr2P = v;
	
	}
	fclose(pFile);

	for(int i=0; i<3; i++){	
		if(calPhys.boolRunDepPar[i]){
			setRunDepPar(calPhys.fileRunDepPar[i], &beam[i], &tar, &lej, &hej);// setting run dependent parameters.
			EBAC = beam[i].EBAC;
			MBeam = beam[i].mass;
			kBF = MFoil/MBeam;
			mA = MBeam; //Beam mass //Reassigned in HandleBOR
			ma = tar.mass;
			mb = lej.mass; //Light ejectile mass
			mB = hej.mass;
	
			printf("Beam energy: %f\n", beam[i].energy);
			printf("Target thickness: %f\n",tar.thickness);
		
			if (grIH2[i]){
				beam[i].energy = beam[i].energy-eloss(beam[i].energy,tar.thickness/2.,grIH2[i]);  
				printf("Energy loss in half target: %f\n" ,eloss(beam[i].energy,tar.thickness/2.,grIH2[i]));
			}
			else printf("Energy loss in target not specified. EBeam=EBAC");

			beam[i].momentum = sqrt(beam[i].energy*beam[i].energy+2.*beam[i].energy*beam[i].mass);//beam momentum
			beam[i].beta = beam[i].momentum/(beam[i].energy + beam[i].mass + tar.mass);
			beam[i].gamma = 1./sqrt(1.-beam[i].beta*beam[i].beta);
			
			EBeam= beam[i].energy;
			PA = beam[i].momentum;//beam momentum
			betaCM = beam[i].beta;
			gammaCM = beam[i].gamma;
	
			printf("Resulting energy: %f\n",EBeam);
			printf("Beam momentum: %f\n",PA);
			printf("Beam beta: %f\tBeam gamma: %f\n",betaCM,gammaCM);
			printf("MBeam: %f\t MFoil: %f\t kBF: %f\n",MBeam,MFoil,kBF);
			printf("beam mass: %f\ttarget mass: %f\n",mA,ma);
			printf("heavy ejectile mass: %f\tlight ejectile mass: %f\n",mB,mb);
			
			//a1 = 1.+mb/mB;
			//a2= 1.-mA/mB;
		}
	}
		printf("End of HandleBOR_Physics\n");
} //HandleBOR_Physics

void HandlePHYSICS(IDet *det)
{
 	//  printf("calling HandlePHYSICS %d %d \n", gCSInitems, gPHOSnitems);
  	//TTree *eventPhy = new TTree ("eventPhy","Physics event"); //Physics event
  	//eventPhy->Branch("fLp",&fLp,"fLp/F"); //Light particle detected at YY1 and CsI
  	cosTheta = cos(TMath::DegToRad()* (det->TSdTheta));

	// Selecting an incoming isotope:	
 	if (calPhys.boolICGates==kFALSE) nGate=0;   
 	else if (calPhys.boolICGates==kTRUE) {  
		for(int i=0; i<3; i++){  
 			if ((det->TICEnergy>beam[i].ICmin) && (det->TICEnergy<beam[i].ICmax)) { 
				nGate = i; 
  			}
		}
	} 

	EBAC = beam[nGate].EBAC;
	MBeam = beam[nGate].mass;
	kBF = MFoil/MBeam;
	mA = MBeam; 
	mB = MBeam;
	EBeam= beam[nGate].energy;
	betaCM = beam[nGate].beta;
	gammaCM = beam[nGate].gamma;
	PA = beam[nGate].momentum; // beam momentum
	// betaCM = sqrt(EBeam*EBeam+2.*EBeam*MBeam)/(EBeam + MBeam + mb);
	// gammaCM = 1./sqrt(1.-betaCM*betaCM);
	// PA = sqrt(EBeam*EBeam+2.*EBeam*MBeam);//beam momentum
	
//---------------------------------------------------------
//	else  if (MBeam == M18O)//18O beam
//    {
//      	//std::cout << "MBeam == M18O" << std::endl; 
//		det->TSdETot = det->TSd2rEnergyCal+det->TSd1rEnergy; //temporarly set to calculated value
//
//	}
//
 	//adding dead layer energy losses
 	//Sd2 ring side
  	if (grIB[nGate]) energy = det->TSd2rEnergy+elossFi(det->TSd2rEnergy,0.1*2.35*0.5/cosTheta,grIB[nGate]); //boron junction implant
  	if (grIAl[nGate]) energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,grIAl[nGate]); //first metal
  	if (grISiO2[nGate]) energy = energy+elossFi(energy,0.1*2.65*2.5/cosTheta,grISiO2[nGate]); //SiO2
  	if (grIAl[nGate]) energy = energy+elossFi(energy,0.1*2.7*1.5/cosTheta,grIAl[nGate]); //second metal
   	//Sd1 ring side
  	if (grIAl[nGate]) energy = energy+elossFi(energy,0.1*2.7*1.5/cosTheta,grIAl[nGate]); //second metal
  	if (grISiO2[nGate]) energy = energy+elossFi(energy,0.1*2.65*2.5/cosTheta,grISiO2[nGate]); //SiO2
  	if (grIAl[nGate]) energy = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,grIAl[nGate]); //first metal
   	energy = energy + det->TSd1rEnergy;// energy lost and measured in Sd1
//	if(MBeam == M20Na){
//		det->TSd1rEnergy = det->TSd1rEnergy+elossFi(energy,0.1*1.822*0.5/cosTheta,grIP[nGate]);    // JSR_july_2 // why???
//		det->TSd1rEnergy = det->TSd1rEnergy+elossFi(energy,0.1*2.7*0.3/cosTheta,grIAl[nGate]);      //JSR_JULY_2
//   	}
	//sector side
	if (grIP[nGate]) energy = energy+elossFi(energy,0.1*1.822*0.5/cosTheta,grIP[nGate]); //phosphorus implant
	if (grIAl[nGate]) det->TSdETot = energy+elossFi(energy,0.1*2.7*0.3/cosTheta,grIAl[nGate]); //metal
    
	PResid = sqrt(2.*det->TSdETot*MBeam);     //Beam momentum in MeV/c
	A = kBF-1.;                              //Quadratic equation parameters
    B = 2.0*PResid* cos(TMath::DegToRad()*det->TSdTheta);
    C = -1.*(kBF+1)*PResid*PResid; 
    if (A!=0)
    PBeam = (sqrt(B*B-4.*A*C)-B)/(2*A);
  	//cout<<"reaching here"<<endl;
  	//to calculate residue energy from beam
 
    // PBeam = sqrt(2.*EBAC*MBeam); //Beam momentum in MeV/c
    // A = kBF+1.; //Quadratic equation parameters
    // B = -2.0*PBeam* cos(TMath::DegToRad()*det->TSdTheta);
    // C = -1.*(kBF-1)*PBeam*PBeam; 
    // if (A!=0)
    // PResid = (sqrt(B*B-4.*A*C)-B)/(2*A);

    det->TBE = PBeam*PBeam/(2.*MBeam);
   
    //  det->TSdETot = EBAC; //tk remove
    //det->TSdThetaCM = TMath::RadToDeg()*atan(sin(TMath::DegToRad()*det->TSdTheta)*PResid/(PResid*cos(TMath::DegToRad()*det->TSdTheta)-PBeam*MBeam/(MBeam+MFoil)));

    // printf("thetaCM: %f\n",det->TSdThetaCM);
    if (grIAg[nGate]) det->TBE=  det->TBE + elossFi(det->TSdETot,foilTh/2., grIAg[nGate]); //energy loss from the end of H2 to the center of Ag.

    det->TSdThetaCM = TMath::RadToDeg()*atan(tan(TMath::DegToRad()*det->TSdTheta)/sqrt(gammaCM-gammaCM*betaCM*(MBeam+det->TBE)/(PBeam*cos(TMath::DegToRad()*det->TSdTheta))));// check if this is still correct for H2 target tk
 	if (grIAg[nGate]) det->TBE=  det->TSdETot + elossFi(det->TSdETot,5.2365/2., grIAg[nGate]); //energy loss from the end of H2 to the center of Ag.
  // cout<<"reaching here 2"<<endl;
   
 // Calculate Q-value from YY1 and CsI// 
	if (((deuterons->IsInside(det->TCsI2Energy[0],det->TYdEnergy[0]*cos(det->TYdTheta[0]*0.01745329)))&& (mb == tar.mass)) && ((det->TYdEnergy[0]>0.2)  && (det->TCsI2Energy[0] >0.6 )&& (mb== tar.mass))) {    //check if in the proton/deuteron gate
//	if (((deuterons->IsInside(det->TCsI2Energy,det->TYdEnergy*cos(det->TYdTheta*0.01745329)))&& (mb == M2H)) && ((det->TYdEnergy>0.2)  && (det->TCsI2Energy >0.6 )&& (mb== M2H))) {    //check if in the proton/deuteron gate

    ECsI= det->TCsI2Energy[0];
    if( ECsI < 0.6){ //approx pedestal vaule // Should probably be changed ?
      	ECsI=0;
	}

  	thetaR =( atan((Yd1rP+((det->TYdRing[0]+1)*(Yd2rP-Yd1rP)/16))/YdDistanceP) + atan((Yd1rP+((det->TYdRing[0])*(Yd2rP-Yd1rP)/16))/YdDistanceP) )/2.;
  	thetaD = thetaR*TMath::RadToDeg();

  	if (mb == tar.mass) //proton energy loss in dead layers between YY1 and CsI                                                                                       
  	//if (mb == M2H) //proton energy loss in dead layers between YY1 and CsI                                                                                       
    {
      	if (grdMy)  ECsI= ECsI+elossFi(ECsI,0.1*1.4*6./cos(thetaR),grdMy); //Mylar                                                                                  
      	else std::cout << "grdMy doesn't exist"<< std::endl;
      	if (grdAl)  ECsI= ECsI+elossFi(ECsI,0.1*2.702*0.3/cos(thetaR),grdAl); //0.3 u Al                                                                            
      	else std::cout << "grdAl doesn't exist"<< std::endl;
      	if (grdP)  ECsI= ECsI+elossFi(ECsI,0.1*1.88219*0.1/cos(thetaR),grdP); // 0.1Phosphorus                                                                      
      	else std::cout << "grdP doesn't exist"<< std::endl;
    }

//---------------YY1 Energy-----------------------------------------------------------------
	EYY1 = det->TYdEnergy[0];
  	//  cout<<"EYY ENERGY IS"<<EYY1<<endl;
 	if(useYCalc){
		if (mb == tar.mass){
		//if (mb == M2H){
  			if (grdSi)  Eb= ECsI + elossFi(ECsI,0.1*2.32*YdThickness[det->TYdNo[0]]/cos(thetaR),grdSi); //Energy loss in YY1 detector // Why calculate a value that you have measured ????? MH
			else std::cout << "grdSi doesn't exist"<< std::endl;
 		}
	}//useYCalc
	else Eb= Eb+EYY1; //use measured Yd // change june28

   	if (mb == tar.mass){
   	//if (mb == M2H){
      	if (grdSi)  Eb= Eb+elossFi(Eb,0.1*2.32*0.35/cos(thetaR),grdSi); //0.3 u Al + 1 um B equivalent in 0.35 um Si                                                            
    	else std::cout << "grdSi doesn't exist"<< std::endl;
    	if (grdD2)  Eb= Eb+elossFi(Eb,tar.thickness/2./cos(thetaR),grdD2); //deuteron energy  in mid target midtarget                                                                             
     	else std::cout << "grdD2 doesn't exist"<< std::endl;
	}

	Pb = sqrt(Eb*Eb+2.*Eb*mb);

//	My eqs  
	Double_t Pby = Pb*sin(thetaR);
	Double_t Pbxcm = gammaCM*betaCM*(Eb+mb)- gammaCM*Pb*cos(thetaR);

// 	cout<<"Beam energy used for QValue is "<< EBeam <<endl;
//  cout<<"deuteron energy at interaction at angle "<< thetaD<< " is equal ="<<Eb<<endl;
 	Q = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb)+2.*PA*Pb*cos(thetaR)+2.*(EBeam+mA+ma-Eb-mb)*ma);  //Alisher's equation 

//	EB = EBeam+mA+ma-Eb-mb; // rel. kinetic energy of recoil.
//	PB = sqrt(PA*PA+Pb*Pb-2.*PA*Pb*cos(thetaR));// momentum of the recoil
//	mB = sqrt(EB^2-PB^2);
//	cout << "Q= " << Q <<endl;   
//	Q = mA+ma-mb-sqrt((EBeam+mA+ma-Eb-mb)*(EBeam+mA+ma-Eb-mb)-(PA*PA+Pb*Pb-2.*PA*Pb*cos(thetaR))); //Equivalent to the previous equation
//	cout << "Q= " << Q <<endl;   
  	det->QValue = Q;
  	IrisEvent->fQv = det->QValue;
  	IrisEvent->fPart.fE=Eb;// proton (deuteron) energy after reaction
	IrisEvent->fThetacm = TMath::RadToDeg()*atan(Pby/Pbxcm);
	if (IrisEvent->fThetacm < 0) 
  		IrisEvent->fThetacm = IrisEvent->fThetacm+180.;
 		#ifdef csv_output
			csv_file << det->QValue <<"," << det->TYdEnergy<<","<<IrisEvent->fPart.fThetacm <<","<< det->TYdTheta<< endl;
		#endif
	}//isinside proton (deuteron) gate
  	IrisEvent->fEYY1 = det->TYdEnergy[0];

	tree->Fill();
}
 //HandlePhysics

Double_t Qvalue (Double_t ma, Double_t mA, Double_t mb,Double_t EBeam, Double_t Eb, Double_t thetaR)
{
  Double_t PA = sqrt(EBeam*EBeam+2.*EBeam*mA);
  Double_t Pb = sqrt(Eb*Eb+2.*Eb*mb);
  return mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb)+2.*PA*Pb*cos(thetaR)+2.*(EBeam+mA+ma-Eb-mb)*ma);
}

void HandleEOR_PHYSICS(int run, int time)
{
 #ifdef csv_output
  csv_file.close();
#endif
  printf(" in Physics EOR\n");
}
 
