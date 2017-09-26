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

extern TFile* treeFile;
extern TTree* tree;

#include "VirtualOdb.h"

extern VirtualOdb* gOdb;

// end

