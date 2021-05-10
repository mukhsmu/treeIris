//
// Global variables for the ROOT analyzer
//
// Name: Globals.h
//
// $Id$
//

// Run parameters
extern int  gRunNumber;
extern int  gEventNumber;
extern bool gIsRunning;
extern bool gIsPedestalsRun;
extern bool gIsOffline;

#include <TFile.h>
#include <TTree.h>
extern TFile* treeFile;
extern TTree* tree;

extern bool gRandomise;

//#include "VirtualOdb.h"
//extern VirtualOdb* gOdb;

//VirtualOdb is obsolete and removed from rootana
#include "mvodb.h"
extern MVOdb *gOdb;

// end

