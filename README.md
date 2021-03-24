# treeIris #
Unpacking software for the IRIS experiment @ TRIUMF


## Download and Installation ##
	
treeIris requires both ROOT (https://root.cern.ch), ROOTANA (https://bitbucket.org/tmidas/rootana) and iris-daqtools (https://bitbucket.org/ttriumfdaq/iris-daqtools).

To install treeIris use

	git clone https://github.com/iris-triumf/treeIris

change variable "ROOTANA" in the Makefile to the path of your ROOTANA installation, 
as well as IRISDIR to the directory where iris-daqtools is installed, and
compile treeIris with

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

The resulting output will be a root file containing the TTree Iris. For each detector the energy (TXxEnergy) and, if requested the raw ADC values (TXxADC) are stored. If the detector has more than one readout channel, the channel number is stored as well (TXxChannel).

For the segmented detectors (YY1, S3, CsI) also the multiplicity (TXxMul), i.e. how many segments were hit, is stored. All other variables are in this case vectors with size TXxMul, and are sorted in descending order by the energy of the hit. E.g., requesting TYdEnergy[0] and TYdChannel[0] will give you the energy and detector segment of the YY1 hit with the highest energy for an event. 

Also for segmented detectors, a few variables depending on the segmentation are stored. For the YY1 detectors, the hit YY1 detector and the number of the hit ring is stored. Also, the reconstructed angle is stored for each of the detectors, depending on the segmentation this can be TXxTheta or TXxPhi.

## Detector Names in root-Tree ##
	
Ionization Chamber:			TIC
Upstream YY1:				TYu
Downstream YY1: 			TYd
CsI, 1st readout:			TCsI1
CsI, 2nd readout:			TCsI2
Upstream S3, ring side:		TSur
Upstream S3, sector side:	TSus
Downstream S3, ring side:	TSd1r
Downstream S3, sector side:	TSd1s
Downstream S3, ring side:	TSd2r
Downstream S3, sector side:	TSd2s
Surface Barrier:			TSSB
Scintillator:				TSc
TRIFIC:						TTr (Was only used in 80Ga test beam)
Zero-degree DSSD, x-strips:	TZdx (Was only used in S1766)
Zero-degree DSSD, y-strips:	TZdy (Was only used in S1766)

## Optional conversion of TDC data and time calibration ##

TDC data can be included in the ROOT output file if treeIris is run with the
'-tdc' option. For detector XXX, the branch TXXXTDC holds the raw TDC
measurement for the channel specified in TXXXTChannel. The TDC signal from the
trigger is stored in TICTDC[0] (i.e. the measurement from TICTChannel==1) and
some jitter can be removed by subtracting this value from the raw TDC
measurements from the other channels. Furthermore, a list of time calibration
files can be provided by running treeIris with '-tdc=/path/to/tdc/config.txt'.
This file must have the format
```
PATH=/path/to/folder/
YD=calibrationYd.txt
SUR=calibrationSur.txt
...
...
```
i.e. much the same format as the main configuration file. Instead of including
PATH, the full file paths can also be specified on the individual lines. The
actual calibration files are two-column text files containing the calibration
coefficients, like
```
slope  delay
0.19354  1948
0.19325  1957
...
...
```
The slope indicates ns per TDC bucket and the delay is the delay between the
trigger and the channel in TDC buckets. The nominal slopes, when running with
'200ps' resolution option in the TDCs, are around 0.195ns/bucket.


## Charge-Sharing Correction ##

Charge-sharing means that a particle deposits in energy not in a single strip of a segmented detector but in several neighbouring ones. In order to correct for this, treeIris will check if reported hits in a single event occured in neighbouring strips of the detector, and sum the energies up. This happens after the hits have been sorted by energy. The lower energy hit gets removed from the vector, the channel number of the hit that was removed is stored in the branch TXxNeighbour.
