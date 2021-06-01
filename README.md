# GoddessAnalysis
## Installation
```
git clone https://github.com/rajghim/GoddessAnalysis.git
cd GoddessAnalysis
mkdir Output
mkdir -p build && cd build
cmake ..
make
```

## Running
After installing, configure the calibrations and configuration described in Calibration and Configuration sections. Three executables (./Q5Gamma, ./SX3Gamma, ./Downstream) are created. Run one of the executables to analyze corresponding detector section.
```
./Q5Gamma (This is to analyze upstream QQQ5 detectors)
./SX3Gamma (This is to analyze upstream SX3 detectors)
./Downstream (This is to analyze downstream barrel telescope)
```

## Calibration
Update the calibration files in "Calibrations" directory.
### BB10EnCal.dat
The columns are:

| Gains | Offsets |
| --- | --- |

### UpQQQ5EnCal.dat
The columns are:

| Gains | Offsets |
| --- | --- |

### (Up/Down)SX3Pedestals.dat
The columns are: 

| Left Pedestals | Right Pedestals | Back Pedestals|
| --- | --- | --- |

### (Up/Down)SX3Gains.dat
The column is:

| SX3 Gains |
| --- |

### (Up/Down)SX3EnCal.dat
The columns are:

| SX3 Front Gains | SX3 Front Offsets | SX3 Back Gains | SX3 Back Offsets |
| --- | --- | --- | ---|

### (Up/Down)SX3PosCal.dat
The columns are: 

| SX3 Left PosCal | SX3 Right PosCal |
| --- | --- |

## Configuration
GoddessAnalysis uses config.json as the main configuration file. As of right now, the config.json is in the following form:
```
{
	"InputPath" : "/mnt/e/goddessSort/Output/",
	"InputFilePrefix": "Run",
	"ICCutPath" : "/mnt/e/Analysis/Analysis Code/Cuts/ICCuts/",
	"ICCutPrefix" : "Pcut_P30dp",
	"OutputPath" : "/mnt/d/Nuclear Phy/GoddessAnalysis/Output/",
	"OutputFilePrefix" : "Goddess30P",
	"Runs":["0048","0049","0050","0078","0079","0081","0083","0084","0088","0089",
			"0102","0103","0104","0105"],
	"BeamMass" : 29.97,
	"TargetMass" : 2.014,
	"EjectileMass" : 1.0072,
	"RecoilMass" : 30.9737,
	"GroundStateQValue" : 10.86,
	"BeamEnergy" : 240,
	"GretinaOffsetX1" : 7.061,
	"GretinaOffsetX2" : 12.37,
	"GretinaOffsetY" : 0.13,
	"GretinaOffsetZ" : 0.2325,
	"ClusterAngle" : 15
}
```

- **InputPath**
	```
	The path to directory containing all the runs.
	```
- **InputFilePrefix**
	```
	The file prefix for the given runs. For a file title 'Run0048_combined.root', this should be 'Run'.
	```
- **ICCutPath**	
	```
	The path to directory containing all the IC cut files. Note that GoddessAnalysis requires individual ic cut file for each run.
	```
- **ICCutFilePrefix**
	```
	The cut file prefix for the given runs. For a file title 'Pcut_p30dp0048.root', this should be 'Pcut_p30dp'
	```
- **OutputPath**
	```
	The path to directory to store output .root files.
	```
- **OutputFilePrefix**	
	```
	The prefix for output file containing output histograms.
	```
- **Runs**
	```
	List all the runs you want to analyze. For example: ["0048","0230"]
	```
- **BeamMass, TargetMass, EjectileMass, RecoilMass, GroundStateQValue**
	```
	Masses of Beam, Target, Ejectile and Recoils respectively (in amu) and Ground State Q-Value (in MeV)
	```
- **GretinaOffsets**
	```
	Gretina offsets in X, Y, and Z axes (in mm)
	```
- **ClusterAngle**
	```
	Maximum angle for cluster addback	
	```

## IC Cut Procedure
Example for **Run0081.root**
1. Open the root file:
	```
	root Run0081.root
	```
2. Open the TBrowser:
	```
	TBrowser a
	```
3. Open icdE VS icE:
	```
	data->Draw("icdE:icE>>(1024,0,4096,1024,0,4096)","","colz")
	```	
4. Select "logz" and zoom in
5. Draw the graphical cut by selecting cut tool from view > Toolbar
6. Use TCutG
	```
	TCUtG *mycut
	mycut = (TCutG*)gROOT->GetListofSpecials()->FindObject("CUTG")
	mycut->SetName("Pcut81")
	```
7. Create Cut File and Write
	```
	TFile *fc = new TFile("Pcut_p30dp0073.root","recreate")
	fc->cd()
	mycut->Write()
	fc->Close()
	.q
	```
	
	


