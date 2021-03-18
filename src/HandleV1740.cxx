#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <assert.h>
#include <signal.h>
#include "TMidasEvent.h"
#include <TFolder.h>
#include <TFile.h>
#include "HandleV1740.h"
#include "Globals.h" //tree

#include <cstdint>
//Pulse processing tools from the iris-daqtools library
//#include <cmath>
//#include "TrapezoidShaper.h"
//#include "LeadingEdgeDiscriminator.h"
//#include "AmplitudeAnalyser.h"
//#include "ConstantFractionAnalyser.h"

extern TTree* tree;

using namespace std;

vector<uint16_t> trace;

void HandleV1740(TMidasEvent& event, void* ptr, int nitems)
{
  trace.clear();
  uint32_t *data;
  data = (uint32_t *) ptr;
  //gV1740nitems = nitems;
  
  if (nitems == 0) { printf("V1740 decode; warning empty bank\n"); return;}

  // This decoding is assuming that group mask == 1 and we are just reading out 8 channels.

  // Structure to store waveforms for 8 channels.
  std::vector<std::vector<int> > samples;
  for(int i = 0; i < 8; i++){
    samples.push_back(std::vector<int>());
  }
  
  // Check that size of data bank looks correct;
  int nblocks = (nitems - 4) / 9;
  int remainder = (nitems - 4) % 9;
  if(remainder != 0) printf("V1740 decoder; bank size is unexpected; size of %i does not satisfy (size-4) % 9 == 0\n",nitems);
  
  // save some information from header.
  int eventID = data[2] & 0xfff;
  uint32_t timeTag = data[3];
  
  //printf("V1740_ Evt#:%d timeTag:%x nwords:%d\n", eventID, timeTag, nitems);


  // Loop over all the 9-word blocks
  for(int i = 0; i < nblocks; i++){
    
    int offset = 4 + i*9;
    // write out the decoding for blocks of 9 words explicitly, otherwise hard to understand this bank decoding
    samples[0].push_back(((data[offset] & 0xfff))); 
    samples[0].push_back(((data[offset] & 0xfff000) >> 12 )); 
    samples[0].push_back(((data[offset] & 0xff000000) >> 24) + ((data[offset+1] & 0xf) <<8) ); 

    samples[1].push_back(((data[offset+1] & 0xfff0) >> 4)); 
    samples[1].push_back(((data[offset+1] & 0xfff0000) >> 16 )); 
    samples[1].push_back(((data[offset+1] & 0xf0000000) >> 28) + ((data[offset+2] & 0xff) <<4) ); 

    samples[2].push_back(((data[offset+2] & 0xfff00)>>8)); 
    samples[2].push_back(((data[offset+2] & 0xfff00000) >> 20 )); 

    offset += 3;

    samples[2].push_back(((data[offset] & 0xfff))); 

    samples[3].push_back(((data[offset] & 0xfff000) >> 12 )); 
    samples[3].push_back(((data[offset] & 0xff000000) >> 24) + ((data[offset+1] & 0xf) <<8) ); 
    samples[3].push_back(((data[offset+1] & 0xfff0) >> 4)); 

    samples[4].push_back(((data[offset+1] & 0xfff0000) >> 16 )); 
    samples[4].push_back(((data[offset+1] & 0xf0000000) >> 28) + ((data[offset+2] & 0xff) <<4) ); 
    samples[4].push_back(((data[offset+2] & 0xfff00)>>8)); 

    samples[5].push_back(((data[offset+2] & 0xfff00000) >> 20 )); 

    offset += 3;

    samples[5].push_back(((data[offset] & 0xfff))); 
    samples[5].push_back(((data[offset] & 0xfff000) >> 12 )); 

    samples[6].push_back(((data[offset] & 0xff000000) >> 24) + ((data[offset+1] & 0xf) <<8) ); 
    samples[6].push_back(((data[offset+1] & 0xfff0) >> 4)); 
    samples[6].push_back(((data[offset+1] & 0xfff0000) >> 16 )); 

    samples[7].push_back(((data[offset+1] & 0xf0000000) >> 28) + ((data[offset+2] & 0xff) <<4) ); 
    samples[7].push_back(((data[offset+2] & 0xfff00)>>8)); 
    samples[7].push_back(((data[offset+2] & 0xfff00000) >> 20 )); 

  }
  
  // save channel 0
  //printf("DT5740: Sample size: %d\n",samples[0].size());
  for(int j = 0; j < samples[0].size(); j++){
    trace.push_back(samples[0][j]);
    //V1730Channel0->SetBinContent(j+1,samples[0][j]);
  }
}
//---------------------------------------------------------------------------------
void HandleBOR_V1740(int file)
{
  printf("V1740 BOR\n");	  
  if(file==0) tree->Branch("wfd",&trace);
	else tree->SetBranchAddress("wfd",&trace);
	printf("Finished HandleBOR_V1740\n");
}
//---------------------------------------------------------------------------------
void HandleEOR_V1740()
{
	printf(" in V1740 EOR\n");
}
