#include "CalibPHYSICS.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>

CalibPHYSICS::CalibPHYSICS(){
  //CalibPHYSICS::Class()->IgnoreTObjectStreamer;
  CalibPHYSICS::Clear();
}

void CalibPHYSICS::ReadFilenames(char* line)
{
	bool expect_val=true;
	char *from=line;
	char *to=line;
	while (*from) {
		if (*from>32) {*to=*from;to++;}
		from++;
	}
	*to=0;
	if (*line==0) return; // line is empty
	
	char* val=strchr(line,'=');
	if (!val){
		val=strchr(line, '!');
		expect_val=false;
	}
	if (!val) printf("Missing = or ! in input file, line: '%s'\n",line);
	*val=0;
	
	// trim param name
	char* trm=val-1;
	while (*trm<=32) *(trm--)=0;
	
	val++;
	if (*val==0 && expect_val) printf("Value missing for parameter %s",line);

	char cval[256];	
	TString strval;
	sscanf(val,"%s",cval);
	strval=cval;
	
	//	parameter of type string:
	if (strcmp(line,"PATH")==0){ 
		installPath = strval;
		fileGates = strval;
		fileELoss = strval;
		filePdedx = strval;
		fileDdedx = strval;
		for(int i=0; i<3; i++){
			fileIdedx[i] = strval;
			fileRunDepPar[i] = strval;
		}
		fileGeometry = strval;
	}

	if (strcmp(line,"GATES")==0){
		fileGates += strval;
		boolGates = kTRUE;
	}
	if (strcmp(line,"ELOSS")==0){
		fileELoss += strval;
		boolELoss = kTRUE;
	}
	if (strcmp(line,"DEDX_P")==0){
		filePdedx += strval;
		boolPdedx = kTRUE;
	}
	if (strcmp(line,"DEDX_D")==0){
		fileDdedx += strval;
		boolDdedx = kTRUE;
	}
	if (strcmp(line,"DEDX_I1")==0){
		fileIdedx[0] += strval;
		boolIdedx[0] = kTRUE;
	}
	if (strcmp(line,"DEDX_I2")==0){
		fileIdedx[1] += strval;
		boolIdedx[1] = kTRUE;
	}
	if (strcmp(line,"DEDX_I3")==0){
		fileIdedx[2] += strval;
		boolIdedx[2] = kTRUE;
	}
	if (strcmp(line,"RUNDEPPAR1")==0){
		fileRunDepPar[0] += strval;
		boolRunDepPar[0] = kTRUE;
	}
	if (strcmp(line,"RUNDEPPAR2")==0){
		fileRunDepPar[1] += strval;
		boolRunDepPar[1] = kTRUE;
	}
	if (strcmp(line,"RUNDEPPAR3")==0){
		fileRunDepPar[2] += strval;
		boolRunDepPar[2] = kTRUE;
	}
	if (strcmp(line,"USE_ICGATES")==0){
		boolICGates = kTRUE;
	}
	if (strcmp(line,"GEOMETRY")==0){
		fileGeometry += strval;
		boolGeometry = kTRUE;
	}
}

void CalibPHYSICS::Load(TString filename){
	boolICGates=kFALSE;
	boolGates=kFALSE;
	boolELoss=kFALSE;
	boolPdedx=kFALSE;
	boolRunDepPar[0]=kFALSE;
	boolRunDepPar[1]=kFALSE;
	boolRunDepPar[2]=kFALSE;
	boolGeometry=kFALSE;

	char line[256];
	FILE* file=fopen(filename.Data(),"rb");
	if (!file)
	{
		printf("Cannot open config file '%s' for reading. Stop.\n",filename.Data());
		exit(0);
	}
	
	printf("Reading config file '%s'\n",filename.Data());
	
	while (!feof(file))
	{
		if (!fgets(line,256,file)) break;
		printf("%s",line);
		// skip leading white spaces
		char* ptr=line;
		while ((*ptr>0) && (*ptr<32)) ptr++;
		//printf("%s\n",ptr[0]);
		switch (ptr[0])
		{
			case 0   :
			case '#' :
			case '/' :  continue;
			default  :  ReadFilenames(ptr);
		}
	}
	fclose(file);
	file=NULL;
}

void CalibPHYSICS::Print(){
	printf("\n********************************\n");
	printf("** PHYSICS calibration files: **\n");
	if(boolGeometry)	printf("Experiment geometry: %s\n",fileGeometry.data());
	else	printf("No experiment geometry specified.\n");
	if(boolELoss)	printf("Energy loss root file: %s\n",fileELoss.data());
	else	printf("No energy loss root file specified.\n");
	if(boolPdedx)	printf("Proton energy loss graphs: %s\n",filePdedx.data());
	else	printf("No proton energy loss graphs specified.\n");
	if(boolDdedx)	printf("Deuteron energy loss graphs: %s\n",fileDdedx.data());
	else	printf("No deuteron energy loss graphs specified.\n");
	if(boolIdedx[0])	printf("Ion energy loss graphs: %s\n",fileIdedx[0].data());
	else	printf("No ion energy loss graphs specified.\n");
	if(boolIdedx[1])	printf("Ion energy loss graphs: %s\n",fileIdedx[1].data());
	else	printf("No ion energy loss graphs specified.\n");
	if(boolIdedx[2])	printf("Ion energy loss graphs: %s\n",fileIdedx[2].data());
	else	printf("No ion energy loss graphs specified.\n");
	if(boolGates)	printf("Deuteron/proton gates: %s\n",fileGates.data());
	else	printf("No deuteron/proton gates specified.\n");
	if(boolRunDepPar[0])	printf("Run dependant parameters: %s\n",fileRunDepPar[0].data());
	else	printf("No run dependant parameters specified.\n");
	if(boolRunDepPar[1])	printf("Run dependant parameters: %s\n",fileRunDepPar[1].data());
	else	printf("No run dependant parameters specified.\n");
	if(boolRunDepPar[2])	printf("Run dependant parameters: %s\n",fileRunDepPar[2].data());
	else	printf("No run dependant parameters specified.\n");
	if(boolICGates)	printf("Using ionization chamber for incoming ID.\n");
	else	printf("No incoming ID using the ionization chamber.\n");
	printf("********************************\n\n");
}

void CalibPHYSICS::Clear(){
}
