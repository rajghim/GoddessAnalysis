#define Analysis_cxx
#include "Analysis.h"
#include "EnergyLoss.h"
#include <TRandom3.h>
#include <string>
#include "EffectiveThickness.h"
#include "Utilities.h"
#include "json/json.h"

//MAIN ANALYSIS LOOP
void Analysis::Loop() {
	DownstreamWelcomeMessage();
	
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    
    //Read and Parse config.json
    Json::Value config;
    std::ifstream config_stream("config.json");
    ASSERT_WITH_MESSAGE(config_stream.is_open(),
                        "Could not find 'config.json'\n");
    config_stream >> config;
    config_stream.close();

    TString ICCutPath = config["ICCutPath"].asString();
    TString ICCutPrefix = config["ICCutPrefix"].asString();
	Float_t Ma = config["BeamMass"].asFloat();
	Float_t Mx = config["TargetMass"].asFloat();
	Float_t Mb = config["EjectileMass"].asFloat();
	Float_t My = config["RecoilMass"].asFloat();
	Float_t Qgs = config["GroundStateQValue"].asFloat();
	Float_t Ta = config["BeamEnergy"].asFloat();
    Float_t GretinaOffsetX1 = config["GretinaOffsetX1"].asFloat(); //Beam Left
    Float_t GretinaOffsetX2 = config["GretinaOffsetX2"].asFloat(); //Beam Right
    Float_t GretinaOffsetY = config["GretinaOffsetY"].asFloat(); //ORRUBA 0.13mm higher than GRETINA
    Float_t GretinaOffsetZ = config["GretinaOffsetZ"].asFloat(); //ORRUBA 0.2325mm downstream than GRETINA
    Float_t ClusterAngle = config["ClusterAngle"].asFloat();
    TString OutputPath = config["OutputPath"].asString();
    TString OutputFilePrefix = config["OutputFilePrefix"].asString();

    //Creating the Output File to Store Histogram
    TFile* outputFile = new TFile(OutputPath + OutputFilePrefix + "_Downstream.root", "recreate");

    //OPEN ALL THE INPUT FILES
	//Open the pedestals File (Keep Left Pedestals and Right Pedestals and Back Pedestals in the same file) 
	std::ifstream file;
	file.open("Calibration/DownSX3pedestals.dat");
	Double_t LeftPed[192] = {0};
	Double_t RightPed[192] = {0};
	Double_t BackPed[192] = {0};
	for (Int_t i = 0; i<192;i++){
		file >> LeftPed[i] >> RightPed[i] >> BackPed[i];	
		//std::cout << LeftPed[1] << '\t' << RightPed[5] << '\t' << BackPed[50] << std::endl;
	}
 
	//Open the gains file
	std::ifstream gainfile;
	gainfile.open("Calibration/DownSX3gains.dat");
	Double_t Gains[48] = {0}; 
	for (Int_t i=0; i<48; i++){
			gainfile >>Gains[i];
			//std::cout << Gains[2] << std::endl;
	}

	//Open the Energy Calibration file for Downstream SX3 (Keep both Front and Back Encal Parameters in same file)
	std::ifstream SX3EnCalfile;
	SX3EnCalfile.open("Calibration/DownSX3EnCal.dat");
	Double_t SX3EnCalSlope[192] = {0};
	Double_t SX3EnCalIntercept[192] = {0};
	Double_t SX3SectorEnCalSlope[192] = {0};
	Double_t SX3SectorEnCalIntercept[192] = {0}; 
	for (Int_t i=0; i<192; i++){
			SX3EnCalfile >> SX3EnCalSlope[i] >> SX3EnCalIntercept[i] >> SX3SectorEnCalSlope[i] >> SX3SectorEnCalIntercept[i];
	}

	//Open the Position Calibration file for Downstream SX3
	std::ifstream SX3PosCalfile;
	SX3PosCalfile.open("Calibration/DownSX3PosCal.dat");
	Double_t xone[48] = {0};
	Double_t xtwo[48] = {0}; 
	for (Int_t i=0; i<48; i++){
			SX3PosCalfile >> xone[i] >> xtwo[i];
	}

	//Open the Energy Calibration file for Downstream BB10
	std::ifstream BB10EnCalfile;
	BB10EnCalfile.open("Calibration/BB10EnCal.dat");
	Double_t BB10EnCalSlope[64] = {0};
	Double_t BB10EnCalIntercept[64] = {0};
	for (Int_t i=0; i<64; i++){
		BB10EnCalfile >> BB10EnCalSlope[i] >> BB10EnCalIntercept[i];
	} 	

	//Define Histograms here
    //IC Stuffs
    TH2F* hICdEE = new TH2F("icdEE", "Ionization Chamber dE vs E; E; dE", 500, 0, 4000, 500, 0, 4000);
    TH2F* hICdEE_cut = new TH2F("icdEE_cut", "Ionization Chamber dE vs E (with IC Cut); E; dE", 500, 0, 4000, 500, 0, 4000);	

    //Total Energy histograms
	TH2F *TotEvA = new TH2F("TotEvA", "Total Energy vs Angle",1000,0,180,2000,0,20000);
	TH2F *TotEvA_TDC = new TH2F("TotEvA_TDC","Total Energy vs Angle with TDC",1000,0,180,2000,0,20000);

	//Test Histograms
	TH2F *BB10EvA_TDC= new TH2F("BB10EvA_TDC","BB10 Energy vs Angle with TDC",1000,0,180,2000,0,20000);
	TH2F *SX3EvA_TDC = new TH2F("SX3EvA_TDC","SX3 Energy vs Angle with TDC",1000,0,180,2000,0,20000);
    TH2F *PID1 = new TH2F("PID1","BB10 Energy vs SX3 Energy",2000,0,20000,2000,0,20000);
    TH2F *BB10VSX3 = new TH2F("BB10VSX3","BB10 Det vs SX3 Det",12,0,12,12,0,12);
    TH2F *BB10stripVSX3strip = new TH2F("BB10stripVSX3strip","BB10 strip vs SX3 Strip",10,0,10,10,0,10);
	
    //ENERGY LOSS DEFINITIONS (DEADLAYER AND TARGET)
	EnergyLoss* DeadLayer = new EnergyLoss("ProtonInSi.dat");
	EnergyLoss* Target = new EnergyLoss("ProtonInC2D4.dat");

    //EVENT BY EVENT ANALYSIS LOOP
    int prevRunNumber = -1;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
	
        // Get run number for TCutG and whatever else
        std::string fileName = fChain->GetCurrentFile()->GetName();
        fileName.erase(fileName.begin(), fileName.end() - 18);
        fileName.erase(fileName.end() - 14, fileName.end());
        TString runNumberStr = fileName;
        int runNumber = std::stoi(fileName);

        if(runNumber != prevRunNumber) {
            TFile *cutFile;
            TString cutFileName = ICCutPath + ICCutPrefix + runNumberStr + ".root";
            std::cout << cyan << "Analyzing Run Number: " << magenta << runNumberStr << reset << std::endl;
            cutFile = TFile::Open(cutFileName);
            TString cutName = Form("Pcut%d", runNumber);
            icdEECut = static_cast<TCutG*>(cutFile->Get(cutName));
            prevRunNumber = runNumber;
            cutFile->Close();
        }

        hICdEE->Fill(icE, icdE);
        //Applying the IC cut here
        if(!icdEECut->IsInside(icE, icdE)) continue;
        // Everything after this has the ic cut on dE vs E
	    hICdEE_cut->Fill(icE, icdE);
		

		//Loop over the BB10 Multiplicity
        for(Int_t k=0; k<BB10Mul; k++){

            //BB10 Calibration
            Float_t BB10Energy = (BB10ADC[k] * (BB10EnCalSlope[(BB10Det[k]*8) + BB10Strip[k]] )) + (BB10EnCalIntercept[(BB10Det[k]*8) + BB10Strip[k]]);
            //Loop over the SX3 Multiplicity
            for(Int_t j=0; j<SX3Mul; j++){    
                if ( (SX3Det[j]==2 && BB10Det[k]==1) || (SX3Det[j]==3 && BB10Det[k]==2) || (SX3Det[j]==8 && BB10Det[k]==5) || (SX3Det[j]==9 && BB10Det[k]==6) || (SX3Det[j]==10 && BB10Det[k]==7) ){
                    if ( (SX3Strip[j]==0 && (BB10Strip[k]==0 || BB10Strip[k]==1 || BB10Strip[k]==2)) || (SX3Strip[j]==1 && (BB10Strip[k]==1 || BB10Strip[k]==2 || BB10Strip[k]==3 || BB10Strip[k]==4)) || (SX3Strip[j]==2 && (BB10Strip[k]==3 || BB10Strip[k]==4 || BB10Strip[k]==5 || BB10Strip[k]==6)) || (SX3Strip[j]==3 && (BB10Strip[k]==5 || BB10Strip[k]==6 || BB10Strip[k]==7))){    
			
			            //Front and Back side Pedestals Substracted
			            Float_t SX3RawStripRight = SX3StripRightADC[j] - RightPed[(SX3Det[j]*16)+(SX3Strip[j]*4)+(SX3Sector[j])];
			            Float_t SX3RawStripLeft = SX3StripLeftADC[j] - LeftPed[(SX3Det[j]*16)+(SX3Strip[j]*4)+(SX3Sector[j])];
			            Float_t SX3RawSector = SX3SectorADC[j] - BackPed[(SX3Det[j]*16)+(SX3Strip[j]*4)+(SX3Sector[j])];
			
			            //Front side Gains applied
			            Float_t RawStripLeft = SX3RawStripLeft; 
			            Float_t RawStripRight = -1. * SX3RawStripRight * (Gains[(SX3Det[j]*4)+(SX3Strip[j])]);  
		
			            //Front sideEnergy and Position Calculations			
			            Float_t RawEnergy = RawStripRight + RawStripLeft; //Gain Matched Energy (No Calibration though)
			            Float_t Energy = RawEnergy * SX3EnCalSlope[(SX3Det[j]*16)+(SX3Strip[j]*4)+(SX3Sector[j])] + SX3EnCalIntercept[(SX3Det[j]*16)+(SX3Strip[j]*4)+(SX3Sector[j])]; //Energy Calibrated
			            if (Energy < 0) continue;
			            Float_t RawPosition = ((RawStripRight - RawStripLeft) / Energy ); // Gain matched Position( Energy has been Calibrated but No Calibration in Position)
			            Float_t Position =((RawPosition - xone[(SX3Det[j]*4)+(SX3Strip[j])]) / (xtwo[(SX3Det[j]*4)+(SX3Strip[j])]- xone[(SX3Det[j]*4)+(SX3Strip[j])])) *75.; //Position Calibration applied

			            //Back side Calibration
			            Float_t SectorEnergy = (SX3RawSector * SX3SectorEnCalSlope[(SX3Det[j]*16)+(SX3Strip[j]*4)+(SX3Sector[j])]) + SX3SectorEnCalIntercept[(SX3Det[j]*16)+(SX3Strip[j]*4)+(SX3Sector[j])];

			            // (x,y,z) for all SX3 detectors (proton hit position)
			            Float_t radius = 100.2;// Radial distance to all detectors from center of the barrel = 100.2mm
                        Float_t Sradius[4]={101.3165,100.3246,100.3246,101.3165};//Radial distances for Strip#-1-2-3-4
			            Float_t ProtonX = radius * sin (((SX3Det[j]*30. + ((SX3Strip[j]-1)*atan2(5.,radius)*180./M_PI) + ((SX3Strip[j]-2)*atan2(5.,radius)*180./M_PI)))*M_PI/180.); //5mm=half strip width
			            Float_t ProtonY = radius * cos(((SX3Det[j]*30. + ((SX3Strip[j]-1)*atan2(5.,radius)*180./M_PI) + ((SX3Strip[j]-2)*atan2(5.,radius)*180./M_PI)))*M_PI/180.);
			            Float_t ProtonZ = Position + 5.; //This number should be 3.3
			
			            //Scaling and Shift Factors
			            Float_t ScaleX = 0; // Don't know the value of this right now
			            Float_t ScaleY = 0; //Don't know the value of this right now
			            Float_t ShiftX = 0;
			            Float_t ShiftY = 0;
			            
                        Float_t TargetX = (icPositionX*ScaleX) + ShiftX;
			            Float_t TargetY = (icPositionY*ScaleY) + ShiftY;
			            Float_t TargetZ = TargetY*tan(27.*M_PI/180.);

			            //Find the distance between the proton hit spot and beam spot and find the reaction angle
			            Float_t d = sqrt((pow((ProtonX-TargetX),2.))+(pow((ProtonY-TargetY),2.))); //Distance between proton spot and beam spot
			            //Float_t SX3Angle = (atan2(d,(ProtonZ + IChamZ))*(180./M_PI));
			            Float_t SX3Angle = (atan2(Sradius[SX3Strip[j]],ProtonZ)*(180./M_PI)); //This is the angle without IC

			            Float_t protonEnergy = BB10Energy + SectorEnergy;
                        Double_t protonEnergyTarget = protonEnergy; //Just doing this for now (no energy loss for now)
			            //Double_t protonEnergyTarget = 1000.*( Target->AddBack(protonEnergy/1000., 0.00292097835/fabs(cos((SX3Angle-27.)*M_PI/180.))) );			
			
			            //Fill your histograms here
                        if (!SX3Upstream[j] && protonEnergy>2600.)TotEvA->Fill(SX3Angle,protonEnergyTarget);
                        if (!SX3Upstream[j] && protonEnergy>2600. && tdcIC>700. && tdcIC<1000.)TotEvA_TDC->Fill(SX3Angle,protonEnergyTarget);
			            if (!SX3Upstream[j] && tdcIC>700. && tdcIC<1000.)BB10EvA_TDC->Fill(SX3Angle,BB10Energy);
			            if (!SX3Upstream[j] && tdcIC>700. && tdcIC<1000.)SX3EvA_TDC->Fill(SX3Angle,SectorEnergy);
                        if (!SX3Upstream[j] && tdcIC>700. && tdcIC<1000.)PID1->Fill(SectorEnergy,BB10Energy);
                        if (!SX3Upstream[j] && tdcIC>700. && tdcIC<1000.)BB10VSX3->Fill(BB10Det[k],SX3Det[j]);
                        if (!SX3Upstream[j] && tdcIC>700. && tdcIC<1000.)BB10stripVSX3strip->Fill(BB10Strip[k],SX3Strip[j]);

                    }//End of if statement for Strip arrangement in Telescope
		        }//End of if statement for Detector arrangement in Telescope	
            }//End of Loop Over SX3 Miltiplicity  
        }//End of Loop Over BB10 Miltiplicity      
    }// End of event by event analysis
  
    //WRITING HISTOGRAMS HERE
    outputFile->cd();

    //IC Stuffs
    hICdEE->Write();
    hICdEE_cut->Write();
    
    //Downstream Stuffs
    TotEvA->Write();
    TotEvA_TDC->Write();
    BB10EvA_TDC->Write();
    SX3EvA_TDC->Write();
    PID1->Write();
    BB10VSX3->Write();
    BB10stripVSX3strip->Write();

    outputFile->Close();
}// END OF THE ANALYSIS LOOP 
