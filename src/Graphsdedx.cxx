#include "Graphsdedx.h"
#include <stdio.h>
#include <stdlib.h>
#include "TString.h"

Graphsdedx::Graphsdedx(){
  //Graphsdedx::Class()->IgnoreTObjectStreamer;
  Graphsdedx::Clear();
}

void Graphsdedx::ReadGraphnames(char* line)
{
	char *from=line;
	char *to=line;
	while (*from) {
		if (*from>32) {*to=*from;to++;}
		from++;
	}
	*to=0;
	if (*line==0) return; // line is empty
	
	char* val=strchr(line,'=');
	if (!val) printf("Missing = in input file, line: '%s'\n",line);
	*val=0;
	
	// trim param name
	char* trm=val-1;
	while (*trm<=32) *(trm--)=0;
	
	val++;
	if (*val==0) printf("Value missing for parameter %s",line);

	char cval[256];	
	TString strval;
	sscanf(val,"%s",cval);
	strval=cval;
	
	//	parameter of type string:
	if (strcmp(line,"Ag")==0){
		Ag = strval;
		boolAg = kTRUE;
	}
	if (strcmp(line,"Al")==0){
		Al = strval;
		boolAl = kTRUE;
	}
	if (strcmp(line,"B")==0){
		B = strval;
		boolB = kTRUE;
	}
	if (strcmp(line,"D2")==0){
		D2 = strval;
		boolP = kTRUE;
	}
	if (strcmp(line,"H2")==0){
		H2 = strval;
		boolH2 = kTRUE;
	}
	if (strcmp(line,"My")==0){
		My = strval;
		boolMy = kTRUE;
	}
	if (strcmp(line,"P")==0){
		P = strval;
		boolP = kTRUE;
	} 
	if (strcmp(line,"Si")==0){
		Si = strval;
		boolSi = kTRUE;
	}
	if (strcmp(line,"SiO2")==0){
		SiO2 = strval;
		boolSiO2 = kTRUE;
	} 
}

void Graphsdedx::Load(TString filename){
	boolAg=kFALSE;
	boolAl=kFALSE;
	boolB=kFALSE;
	boolD2=kFALSE;
	boolH2=kFALSE;
	boolMy=kFALSE;
	boolP=kFALSE;
	boolSi=kFALSE;
	boolSiO2=kFALSE;

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
			default  :  ReadGraphnames(ptr);
		}
	}
	fclose(file);
	file=NULL;
}

void Graphsdedx::Print(){
	if(boolAg==kTRUE)	printf("Ag graph: %s\n",Ag.Data());
	else	printf("Ag graph not specified.\n");
	if(boolAl==kTRUE)	printf("Al graph: %s\n",Al.Data());
	else	printf("Al graph not specified.\n");
	if(boolB==kTRUE)	printf("B graph: %s\n",B.Data());
	else	printf("B graph not specified.\n");
	if(boolD2==kTRUE)	printf("D2 graph: %s",P.Data());
	else	printf("D2 graph not specified.\n");
	if(boolH2==kTRUE)	printf("H2 graph: %s\n",H2.Data());
	else	printf("H2 graph not specified.\n");
	if(boolMy==kTRUE)	printf("My graph: %s\n",My.Data());
	else	printf("My graph not specified.\n");
	if(boolP==kTRUE)	printf("P graph: %s\n",P.Data());
	else	printf("P graph not specified.\n");
	if(boolSi==kTRUE)	printf("Si graph: %s\n",Si.Data());
	else	printf("Si graph not specified.\n");
	if(boolSiO2==kTRUE)	printf("SiO2 graph: %s\n",SiO2.Data());
	else	printf("SiO2 graph not specified.\n");
}

void Graphsdedx::Clear(){
}
