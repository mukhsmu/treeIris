// run dependent parameters
#ifndef runDep_H
#define runDep_H
#include <TObject.h>
#include <TClass.h>
#include <string>

class runDep : public TObject {
	public:
		runDep(); 
		virtual ~runDep() {} //! 
	 	double EBAC; //beam energy from accelerator
	 	double energy; //beam energy behind target (also initialized with value from accelerator)
	  	double momentum;// momentum in MeV. Initial value calculated using EBAC, has to be properly set in HandleBOR
	  	double beta;// beam beta in center-of-mass system. Initial value calculated using EBAC, has to be properly set in HandleBOR
	  	double gamma;// beam gamma in center-of-mass system. Initial value calculated using EBAC, has to be properly set in HandleBOR 
		float ICmin; // Lower limit of ionization chamber gate for incoming beam.
		float ICmax; // Upper limit of ionization chamber gate for incoming beam.
		std::string nA;
		std::string na;
		std::string nB;
		std::string nb;

		void setRunDepPar(std::string);
		void Print();
};
#endif
// end
