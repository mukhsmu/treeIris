# treeIris #
Analysis software for the IRIS experiment @ TRIUMF


## Download and Installation ##
	
treeIris requires both ROOT (https://root.cern.ch) and ROOTANA (https://bitbucket.org/tmidas/rootana).

To install treeIris use

	git clone https://github.com/iris-triumf/treeIris

change variable "ROOTANA" in the Makefile to the path of your ROOTANA installation. And compile treeIris with

	make


## Running treeIris ##

From the main treeIris folder, you can run treeIris with
 
	./bin/treeIris -c=/path/to/your/config.file -o=/path/to/to/your/output.file /path/to/your/input.file

The input file has to be a MIDAS file (eg. *.mid or *.mid.gz), the output file will be a root file. Uncalibrated ADC data can be requested with the option '-raw', TDC data with the option '-tdc'.


## The Config File ##

In the main config file the paths to further files which in turn contain detector calibration parameters are defined. The structure is always

	FILETYPE=/path/to/calibration.file

The following files can be specified

	GEOMETRY : Contains information on detector distances and thicknesses, and the target thickness and configuration. 
	SD1R : Pedestal and gain for the ring-side of the first downstream S3 detector.
	SD1S : Pedestal and gain for the sector-side of the first downstream S3 detector.
	SD2R : Pedestal and gain for the ring-side of the second downstream S3 detector.
	SD2S : Pedestal and gain for the sector-side of the second downstream S3 detector.
	SUR : Pedestal and gain for the ring-side of the upstream S3 detector.
	SUS : Pedestal and gain for the ring-side of the upstream S3 detector.
	CSI1 : Pedestal and gain for the first readout of the CsI detector.
	CSI2 : Pedestal and gain for the second  readout of the CsI detector.
	YD : Pedestal and gain for the downstream YY1 detector.
	YU : Pedestal and gain for the upstream YY1 detector.
	TCORR_IC : Defines run-dependent corrections to the energy loss in the IC. 
	TCORR_SI : Defines run-dependent corrections to the S3 calibrations. 

Without any configuration file, the program will generate a file containing uncalibrated data.

## The Output File Structure ##

The resulting output will be a root file containing the TTree Iris. This tree contains variables for the hit multiplicity, the energy, the hit segment, the measured angle, and, if requested, the time, for each detector.
