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
#include <cmath>
//#include "TrapezoidShaper.h"
//#include "LeadingEdgeDiscriminator.h"
//#include "AmplitudeAnalyser.h"
//#include "ConstantFractionAnalyser.h"

extern TTree* tree;

using namespace std;

vector<vector<uint16_t>> traces(16,vector<uint16_t>());
//vector<double> amplitudes;
//vector<double> times;

void HandleV1740(TMidasEvent& event, void* ptr, int nitems)
{
  //Clear old values
  for(unsigned int i=0; i<traces.size(); i++) traces.at(i).clear();
  //trace.clear();
  //amplitudes.clear();
  //times.clear();

  uint32_t *data;
  data = (uint32_t *) ptr;
  //gV1740nitems = nitems;
  
  if (nitems == 0) { printf("V1740 decode; warning empty bank\n"); return;}

  // This decoding is assuming that group mask == 1 and we are just reading out 8 channels.
  // Update: No assumption about group mask anymore.

  // Structure to store waveforms for 8 channels.
  //std::vector<std::vector<int> > samples;
  //for(int i = 0; i < 8; i++){
  //  samples.push_back(std::vector<int>());
  //}
  
  // Check that size of data bank looks correct;
  int nblocks = (nitems - 4) / 9;
  int remainder = (nitems - 4) % 9;
  if(remainder != 0) printf("V1740 decoder; bank size is unexpected; size of %i does not satisfy (size-4) % 9 == 0\n",nitems);
  
  // save some information from header.
  uint8_t groupMask = data[1] & 0xf;
  vector<uint8_t> groups;
  //Which groups participate?
  for(uint8_t i=0; i<4; i++) if((groupMask >> i) & 0x01) groups.push_back(i);
  uint8_t nGroups = groups.size();

  vector<vector<uint16_t>> samples (32,vector<uint16_t>());

  int eventID = data[2] & 0xfff;
  uint32_t timeTag = data[3];
  
  //printf("V1740_ Evt#:%d timeTag:%x nwords:%d\n", eventID, timeTag, nitems);

  for(unsigned int j=0; j<nGroups; j++){
    uint8_t group_j = groups.at(j);
    // Loop over all the 9-word blocks
    for(unsigned int i = j*nblocks/nGroups; i < (j+1)*nblocks/nGroups; i++){
    
      int offset = 4 + i*9;
      // write out the decoding for blocks of 9 words explicitly, otherwise hard to understand this bank decoding
      
      traces[group_j*4+0].push_back(((data[offset] & 0xfff))); 
      traces[group_j*4+0].push_back(((data[offset] & 0xfff000) >> 12 )); 
      traces[group_j*4+0].push_back(((data[offset] & 0xff000000) >> 24) + ((data[offset+1] & 0xf) <<8) ); 
      
      traces[group_j*4+1].push_back(((data[offset+2] & 0xfff00)>>8)); 
      traces[group_j*4+1].push_back(((data[offset+2] & 0xfff00000) >> 20 )); 
      offset += 3;
      traces[group_j*4+1].push_back(((data[offset] & 0xfff)));
      
      traces[group_j*4+2].push_back(((data[offset+1] & 0xfff0000) >> 16 )); 
      traces[group_j*4+2].push_back(((data[offset+1] & 0xf0000000) >> 28) + ((data[offset+2] & 0xff) <<4) ); 
      traces[group_j*4+2].push_back(((data[offset+2] & 0xfff00)>>8));
      
      offset += 3;
          
      traces[group_j*4+3].push_back(((data[offset] & 0xff000000) >> 24) + ((data[offset+1] & 0xf) <<8) ); 
      traces[group_j*4+3].push_back(((data[offset+1] & 0xfff0) >> 4)); 
      traces[group_j*4+3].push_back(((data[offset+1] & 0xfff0000) >> 16 )); 
      
      //samples[group_j*8+0].push_back(((data[offset] & 0xfff))); 
      //samples[group_j*8+0].push_back(((data[offset] & 0xfff000) >> 12 )); 
      //samples[group_j*8+0].push_back(((data[offset] & 0xff000000) >> 24) + ((data[offset+1] & 0xf) <<8) ); 

      //samples[group_j*8+1].push_back(((data[offset+1] & 0xfff0) >> 4)); 
      //samples[group_j*8+1].push_back(((data[offset+1] & 0xfff0000) >> 16 )); 
      //samples[group_j*8+1].push_back(((data[offset+1] & 0xf0000000) >> 28) + ((data[offset+2] & 0xff) <<4) ); 
      
      //samples[group_j*8+2].push_back(((data[offset+2] & 0xfff00)>>8)); 
      //samples[group_j*8+2].push_back(((data[offset+2] & 0xfff00000) >> 20 )); 
      //offset += 3;
      //samples[group_j*8+2].push_back(((data[offset] & 0xfff))); 

      //samples[group_j*8+3].push_back(((data[offset] & 0xfff000) >> 12 )); 
      //samples[group_j*8+3].push_back(((data[offset] & 0xff000000) >> 24) + ((data[offset+1] & 0xf) <<8) ); 
      //samples[group_j*8+3].push_back(((data[offset+1] & 0xfff0) >> 4)); 

      //samples[group_j*8+4].push_back(((data[offset+1] & 0xfff0000) >> 16 )); 
      //samples[group_j*8+4].push_back(((data[offset+1] & 0xf0000000) >> 28) + ((data[offset+2] & 0xff) <<4) ); 
      //samples[group_j*8+4].push_back(((data[offset+2] & 0xfff00)>>8)); 

      //samples[group_j*8+5].push_back(((data[offset+2] & 0xfff00000) >> 20 )); 

      //offset += 3;
      //samples[group_j*8+5].push_back(((data[offset] & 0xfff))); 
      //samples[group_j*8+5].push_back(((data[offset] & 0xfff000) >> 12 )); 

      //samples[group_j*8+6].push_back(((data[offset] & 0xff000000) >> 24) + ((data[offset+1] & 0xf) <<8) ); 
      //samples[group_j*8+6].push_back(((data[offset+1] & 0xfff0) >> 4)); 
      //samples[group_j*8+6].push_back(((data[offset+1] & 0xfff0000) >> 16 )); 

      //samples[group_j*8+7].push_back(((data[offset+1] & 0xf0000000) >> 28) + ((data[offset+2] & 0xff) <<4) ); 
      //samples[group_j*8+7].push_back(((data[offset+2] & 0xfff00)>>8)); 
      //samples[group_j*8+7].push_back(((data[offset+2] & 0xfff00000) >> 20 )); 
    }
  }
  
  //Save data.
  //for(unsigned int i=0; i<samples.size(); i+=2){
  //  for(unsigned int j=0; j<samples[i].size(); j++){
  //    traces.at(i/2).push_back(samples[i][j]);
  //  }
  //}
    
  // save channel 0
  //printf("DT5740: Sample size: %d\n",samples[0].size());
  //for(int j = 0; j < samples[0].size(); j++){
    //trace.push_back(samples[0][j]);
    //V1730Channel0->SetBinContent(j+1,samples[0][j]);
  //}

  //XXX: Analysis of pulse heights and timings (cfg. parameters should come
  //     from a configuration file and not be hard-coded).
  //Shaping parameters:
  /*
  int M = 937; //Corresponds to 15us decay time.
  int k = 12;  //Rise time
  int m = 6;   //Top width
  int N = samples[0].size(); //Length of recorded trace.
  //CFD parameters:
  double att_factor = 0.3;
  int delay = 10;
  //Edge detection parameters:
  int gate = 2*k+m;        //Analysis gate (= dead time).
  double threshold = 50; //Triggering threshold.
  bool negative = false;   //Polarity (false <-> positive polarity).

  //Shape the signal
  TrapezoidShaper shaper(M,k,m);
  double *shaped = shaper.Shape(samples[0].data(),N);  

  //Identify edges in the shaped signal.
  LeadingEdgeDiscriminator discriminator(threshold,negative);
  discriminator.SetDeadTime(gate);
  int hits = discriminator.Analyse(shaped,N);
  std::vector<double> & triggers =  discriminator.GetTriggers();
  
  //We loop over the hits and fill the spectra.
  for(int j=0; j<hits; j++){
    int start_bin = std::floor(triggers.at(j));
    double *ptr = shaped + start_bin;
    //Extract pulse heights.
    AmplitudeAnalyser e_analyser(negative);
    double E = e_analyser.Analyse(ptr,gate);
    amplitudes.push_back(E);
    //Extract pulse timings.
    ConstantFractionAnalyser t_analyser(att_factor,delay,negative);
    double t = t_analyser.Analyse(ptr,gate);
    times.push_back(start_bin+t);
  }
  */
}
//---------------------------------------------------------------------------------
void HandleBOR_V1740(int file)
{
  printf("V1740 BOR\n");	  
  if(file==0){
    for(int i=0; i<16; i++){
      char name[32];
      sprintf(name,"TWfdTrace%02d",2*i);
      tree->Branch(name,&(traces.at(i)));
    //tree->Branch("TWfdADC",&amplitudes);
    //tree->Branch("TWfdTDC",&times);
    }
  }
  else{
	 for(int i=0; i<16; i++){
      char name[32];
      sprintf(name,"TWfdTrace%02d",2*i);
      tree->SetBranchAddress(name,&(traces.at(i)));
    //tree->SetBranchAddress("TWfdADC",&amplitudes);
    //tree->SetBranchAddress("TWfdTDC",&times);
    }
  }
	printf("Finished HandleBOR_V1740\n");
}
//---------------------------------------------------------------------------------
void HandleEOR_V1740()
{
	printf(" in V1740 EOR\n");
}
