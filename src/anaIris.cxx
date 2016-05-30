// treeIris
// Analysis software for the IRIS experiment @TRIUMF
// based on
// ROOT analyzer
// by K.Olchanski

#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <assert.h>
#include <signal.h>
#include <string>

#include "TMidasEvent.h"
#include "TMidasFile.h"
#include "XmlOdb.h"

#include <TTimer.h>
#include <TFile.h>
#include <TTree.h>
#include "Globals.h"

#include "HandleMesytec.h"
#include "HandleV1190.h"
#include "HandleSTAT.h"
#include "HandleScaler.h"
#include "HandlePHYSICS.h"

// Global Variables
int  gRunNumber = 0;
bool gIsPedestalsRun = false;
int  gEventCutoff = 0;
bool gUseTdc = false;
bool gUseRaw = false;
std::string gOutputFile = "";
std::string gCalibFile = "";

char mesbkname[][5] = {"ICA_", "SD2A", "SD1A", "YDA_","SUA_",""};
char tdcbkname[][5] = {"ICT_", "YDT_", "YUT_", "SUT_" ,"SD2T" ,"SD1T"};
char stabkname[][5] = {"VTST", "VTSA", "VTRT", "VTRA"};
char scalbkname[][5] = {"SCAD", "SCAR", "SCAS",""}; //scalers

TFile* treeFile = NULL;
TTree *tree = NULL;
TEvent *IrisEvent = NULL;
VirtualOdb* gOdb = NULL;
FILE* ASCIIYY1 = NULL;
FILE* ASCIICsI1 = NULL;
FILE* ASCIICsI2 = NULL;
FILE* ASCIIIC = NULL;
FILE* ASCIISd1 = NULL;
FILE* ASCIISd2 = NULL;

IDet detec; // calibrated variables from detectors, to be passed to HandlePhysics
IDet *pdet = &detec;
ITdc timeArray;
ITdc *ptdc = &timeArray;

double GetTimeSec() {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec + 0.000001*tv.tv_usec;
}

//--------------------------------------------------------------
void startRun(int transition,int run,int time)
{
	gRunNumber = run;
	gIsPedestalsRun = gOdb->odbReadBool("/experiment/edit on start/Pedestals run");
	printf("Begin run: %d, pedestal run: %d\n", gRunNumber, gIsPedestalsRun);
		
	if(treeFile!=NULL)
	{
	  	treeFile->Write();
	  	treeFile->Close();
	  	treeFile=NULL;
	}  
	
	char filename[1024];
	
	sprintf(filename, "%s", gOutputFile.data());
	treeFile = new TFile(filename,"RECREATE");
	printf("Writing data to %s.\n",gOutputFile.data());
	
	HandleBOR_Mesytec(run, time, pdet, gCalibFile);
	if(gUseTdc) HandleBOR_V1190(run, time, ptdc);
	HandleBOR_PHYSICS(run,time, pdet, gCalibFile); 
	HandleBOR_Scaler(run,time); 
	HandleBOR_STAT(run, time);
}

void endRun(int transition,int run,int time)
{
	gRunNumber = run;
	
	HandleEOR_Mesytec(run, time);
	if(gUseTdc) HandleEOR_V1190(run, time);
	HandleEOR_PHYSICS(run,time);
	HandleEOR_Scaler(run, time); 
	HandleEOR_STAT(run, time);
 
	if (treeFile)
    {
      	tree->AutoSave();
      	treeFile->Close();
      	treeFile = NULL;
    }

  	printf("End of run %d\n",run);
}

//--------------------------------------------------------------
void HandleMidasEvent(TMidasEvent& event)
{
  	void *ptr;
  	int m;
  	int eventId = event.GetEventId();
  	//uint32_t timeStamp = event.GetTimeStamp();
  	//uint32_t serialN = event.GetSerialNumber();
  
  	// clear values here!
	clearDetectors();
  	if ((eventId == 1)) {
    //
    // Mesytec modules
    pdet->Clear();
  	IrisEvent->Clear();
    m=0;
    while (mesbkname[m][0]) { 
      	int size = event.LocateBank(NULL, mesbkname[m], &ptr);
      	//  if (m==2)  printf("ADC bank %s: first: %d %d\n", mesbkname[m], mesbkname[m][0], size); // loop until ascii code for first character of the bank name is empty
      	if (ptr && size) { 
			HandleMesytec(event, ptr, size, m, pdet, gCalibFile, gUseRaw); 
 		} 
      	m++;
 	}
	if(gUseTdc){
		m=0;
  		while (tdcbkname[m][0]) {
			int size = event.LocateBank(NULL, tdcbkname[m], &ptr);
			// printf("TDC bank %s:%d\n", tdcbkname[m], size);
			if (ptr && size) {
		  		HandleV1190(event, ptr, ptdc, size, m); 
			}
			m++;
  		}
	}
	HandlePHYSICS(pdet);
  	}
 	else if ((eventId == 3)) { // Scaler modules
      	m=0;
    	while (scalbkname[m][0]) { 
      		int size = event.LocateBank(NULL, scalbkname[m], &ptr);
      		//  if (m==2)  printf("ADC bank %s: first: %d %d\n", mesbkname[m], mesbkname[m][0], size); // loop until ascii code for first character of the bank name is empty
      		if (ptr && size) { 
				HandleScaler(event, ptr, size, m); 
      		} 
      		m++;
    	}
  	}

#if 0
  // Stat banks
  int size =  event.LocateBank(NULL, stabkname[3], &ptr);
  printf("bank %s:%d\n", tdcbkname[3], size);              
  if (ptr && size) {
    HandleSTAT(event, ptr, size, 0);
  }
#endif

}

//--------------------------------------------------------------
int ProcessMidasFile(/*TApplication*app,*/const char*fname)
{
	char dcfname[100] = "dccp://";
   	strcat(dcfname,fname); // added dccp:// in front of fname, there is seg fault
  	TMidasFile f;
  	bool tryOpen = f.Open(dcfname);
 
  	if (!tryOpen)
    {
    	printf("Cannot open input file \"%s\"\n",dcfname);
      	return -1;
    }

  	int i=0;
  	while (1)
    {
      	TMidasEvent event;
      	if (!f.Read(&event))
		break;

      	int eventId = event.GetEventId();

      	if ((eventId & 0xFFFF) == 0x8000)
		{
	  		// begin run
	  		printf("---- BEGIN RUN ---- \n");
	  		// Load ODB contents from the ODB XML file
	  		if (gOdb) delete gOdb;
	  		gOdb = new XmlOdb(event.GetData(),event.GetDataSize());
	  		startRun(0,event.GetSerialNumber(),0);
		}
      	else if ((eventId & 0xFFFF) == 0x8001)
		{
	  		printf("---- END RUN ---- \n");
	  		// end run
		}
      	else
		{
	  		event.SetBankList();
	  		HandleMidasEvent(event);
		}
      
      	if((i%100)==0)
		{
	   		printf("Processing event %d\r",i);
		}
       	i++;
      	if ((gEventCutoff!=0)&&(i>=gEventCutoff))
		{
	  		printf("Reached event %d, exiting loop.\n",i);
	  		break;
		}
    }
  
  	f.Close();
  	endRun(0,gRunNumber,0);

  	return 0;
}

void help()
{
	printf("\nUsage:\n");
	printf("\n./bin/anaIris [-h] [-c][-o] [-eMaxEvents] [file1 file2 ...]\n");
	printf("\n");
	printf("\t-h: print this help message\n");
	printf("\t-o/-o=/--output=: Path of output file.\n");
	printf("\t-c/-c=/--config=: Path of main configuration file.\n");
	printf("\t-e: Number of events to read from input data files\n");
	printf("\t-nt: Don't read TDCs.\n");
	printf("\n");
	exit(1);
}

// Main function call

int main(int argc, char *argv[])
{
	setbuf(stdout,NULL);
	setbuf(stderr,NULL);
	
	signal(SIGILL,  SIG_DFL);
	signal(SIGBUS,  SIG_DFL);
	signal(SIGSEGV, SIG_DFL);
	
	std::vector<std::string> args;
	for (int i=0; i<argc; i++)
	{
		if (strcmp(argv[i],"-h")==0) help(); // does not return
	    args.push_back(argv[i]);
	}
	
	bool have_output= false;
	bool have_calib= false;

   	for (unsigned int i=1; i<args.size(); i++) // loop over the commandline options
 	{
       	const char* arg = args[i].c_str();
	   
        if (strncmp(arg,"-e",2)==0)  // Event cutoff flag (only applicable in offline mode)
	 		gEventCutoff = atoi(arg+2);
       	else if (strncmp(arg,"-o=",3)==0){  // Output file 
	 		gOutputFile = arg+3;
			have_output=true;
		}
       	else if (strncmp(arg,"-o",2)==0){  // Output file 
	 		gOutputFile = arg+2;
			have_output=true;
		}
      	else if (strncmp(arg,"--output=",9)==0){  // Output file 
	 		gOutputFile = arg+9;
			have_output=true;
		}
		else if (strncmp(arg,"-c=",3)==0){  // Calibration file 
	 		gCalibFile = arg+3;
			have_calib=true;
		}
       	else if (strncmp(arg,"-c",2)==0){  // Calibration file 
	 		gCalibFile = arg+2;
			have_calib=true;
		}
      	else if (strncmp(arg,"--config=",9)==0){  // Calibration file 
	 		gCalibFile = arg+9;
			have_calib=true;
		}
       	else if (strcmp(arg,"-tdc")==0)
	 		gUseTdc=true;
       	else if (strcmp(arg,"-raw")==0)
	 		gUseRaw=true;
       	else if (strcmp(arg,"-h")==0)
	 		help(); // does not return
       	else if (arg[0] == '-')
	 		help(); // does not return
	}
	
	for (unsigned int i=1; i<args.size(); i++)
    {
       	const char* arg = args[i].c_str();

      	if (arg[0] != '-')  
	 	{  
			if(have_output==true&&have_calib==true){	   
				ProcessMidasFile(/*app,*/arg);
			}
			else{
				if(have_output==false){
					printf("No output file specified!\n");
					printf("Use -o=/path/to/your/output-file!\n");
				}
				if(have_calib==false){ 
					printf("No calibration file specified!\n");
					printf("Use -c=/path/to/your/calib-file!\n");
				}
	 		}
     	}
	}
   	return 0;
}

//end
