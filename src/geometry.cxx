//Run dependent parameters in BOR AS
#include "geometry.h"

geometry::geometry(){
}

void geometry::ReadGeometry(std::string filename)
{
	YdDistance = 0.;
	Yd1r = 0.;
	Yd2r = 0.;
	Sd1Distance = 0.;
	Sd2Distance = 0.;
	Sdr1 = 0.;
	Sdr2 = 0.;
	xShift = 0.;
	yShift = 0.;

	char buffer[256];
	FILE * parFile;
	parFile=fopen(filename.data(),"r");

	if (parFile == NULL) {
		perror ("Error opening file");
   	}	
	printf("Reading config file '%s'\n",filename.data());
	
	while (!feof(parFile))
	{
		if (!fgets(buffer,256,parFile)) break;
		printf("%s",buffer);
		
		char* val=strchr(buffer,'=');
		if (!val) printf("Missing = in input parFile, line: '%s'\n",buffer);
		*val=0;	
		
		// trim param name
		char* trm=val-1;
		while (*trm<=32) *(trm--)=0;
	
		val++;
		if (*val==0) printf("Value missing for parameter %s",buffer);
		
		// parse float parameter (if any)
		Double_t v;
		sscanf(val,"%lf",&v);

		if (strcmp(buffer,"YDD")==0)	YdDistance = v;
		if (strcmp(buffer,"YDR1")==0)	Yd1r = v;
		if (strcmp(buffer,"YDR2")==0)	Yd2r = v;
		if (strcmp(buffer,"SD1D")==0)	Sd1Distance = v;
		if (strcmp(buffer,"SD2D")==0)	Sd2Distance = v;
		if (strcmp(buffer,"SDR1")==0)	Sdr1 = v;
		if (strcmp(buffer,"SDR2")==0)	Sdr2 = v;
		if (strcmp(buffer,"XSHIFT")==0)	xShift = v;
		if (strcmp(buffer,"YSHIFT")==0)	yShift = v;
		if (strcmp(buffer,"TTH")==0)	TThickness = v;
	}
	fclose(parFile);
}

void geometry::Print()
{
	printf("\n*** Detector Geometry:\n");		
	printf("Inner radius YY1: %.1f mm\n",Yd1r);
	printf("Outer radius YY1: %.1f mm\n",Yd2r);
	printf("Inner radius S3: %.1f mm\n",Sdr1);
	printf("Outer radius S3: %.1f mm\n",Sdr2);
	printf("Distance YY1: %.1f mm\n",YdDistance);
	printf("Distance 1st S3: %.1f mm\n",Sd1Distance);
	printf("Distance 2nd S3: %.1f mm\n",Sd2Distance);
	printf("Beam shift in x-direction: %.1f mm\n",xShift);
	printf("Beam shift in y-direction: %.1f mm\n",yShift);
	printf("Target thickness: %.2f mg/cm^2\n",TThickness);
}
