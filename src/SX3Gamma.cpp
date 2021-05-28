#define Analysis_cxx
#include "Analysis.h"
#include "EnergyLoss.h"
#include <TRandom3.h>
#include <string>
#include "EffectiveThickness.h"
#include "json/json.h"

//DOPPLER CORRECTION FUNCTION DEFINITION
Float_t dop(Float_t GamAngle, Float_t RecoilBeta, Float_t GamEnergy){
    Float_t RecoilGamma = 1 / (sqrt (1-(RecoilBeta*RecoilBeta) ) );
    Float_t DopplerCorrectedEnergy = GamEnergy * ( 1-RecoilBeta*(cos(GamAngle*M_PI/180.)) ) * RecoilGamma;
    return DopplerCorrectedEnergy;
}

void Analysis::Loop() {
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
    TFile* outputFile = new TFile(OutputPath + OutputFilePrefix + "_SX3.root", "recreate");



	//Open the pedestals File (Keep Left Pedestals and Right Pedestals and Back Pedestals in the same file) 
	std::ifstream file;
	file.open("Calibration/UpSX3pedestals.dat");
	Double_t LeftPed[192] = {0};
	Double_t RightPed[192] = {0};
	Double_t BackPed[192] = {0};
	for (Int_t i = 0; i<192;i++){
		file >> LeftPed[i] >> RightPed[i] >> BackPed[i];	
		//std::cout << LeftPed[1] << '\t' << RightPed[5] << '\t' << BackPed[50] << std::endl;
	}
 
	//Open the gains file
	std::ifstream gainfile;
	gainfile.open("Calibration/UpSX3gains.dat");
	Double_t Gains[48] = {0}; 
	for (Int_t i=0; i<48; i++){
			gainfile >>Gains[i];
			//std::cout << Gains[2] << std::endl;
	}

	//Open the Energy Calibration file for Upstream SX3 (Keep both Front and Back Encal Parameters in same file)
	std::ifstream SX3EnCalfile;
	SX3EnCalfile.open("Calibration/UpSX3EnCal.dat");
	Double_t SX3EnCalSlope[192] = {0};
	Double_t SX3EnCalIntercept[192] = {0};
	Double_t SX3SectorEnCalSlope[192] = {0};
	Double_t SX3SectorEnCalIntercept[192] = {0}; 
	for (Int_t i=0; i<192; i++){
			SX3EnCalfile >> SX3EnCalSlope[i] >> SX3EnCalIntercept[i] >> SX3SectorEnCalSlope[i] >> SX3SectorEnCalIntercept[i];
	}

	//Open the Position Calibration file for Upstream SX3
	std::ifstream SX3PosCalfile;
	SX3PosCalfile.open("Calibration/UpSX3PosCal.dat");
	Double_t xone[48] = {0};
	Double_t xtwo[48] = {0}; 
	for (Int_t i=0; i<48; i++){
			SX3PosCalfile >> xone[i] >> xtwo[i];
	}

	//DEFINE HISTOGRAMS HERE
	//IC stuffs
	TH2F* hICdEE = new TH2F("icdEE", "Ionization Chamber dE vs E; E; dE", 500, 0, 4000, 500, 0, 4000);
    TH2F* hICdEE_cut = new TH2F("icdEE_cut", "Ionization Chamber dE vs E (with IC Cut); E; dE", 500, 0, 4000, 500, 0, 4000);	

	//SX3 Stuffs
	TH2F *SX3EvA = new TH2F("SX3EvA", "SX3 Energy vs Angle",1000,0,180,2000,0,20000);
	TH2F *SX3EvA_TDC = new TH2F("SX3EvA_TDC","SX3 Energy vs Angle with TDC",1000,0,180,2000,0,20000);
	TH1F *SX3Ex = new TH1F("SX3Ex","SX3 Excitation Energy",250,-2,10);
	TH1F *SX3Ex_TDC = new TH1F("SX3Ex_TDC","SX3 Excitation Energy with TDC",250,-2,10);

	//GammaRay Stuffs
	TH2F *GenVQen = new TH2F("GenVQen","Gamma Energy vs Excitation Energy",1000,0,10000,20000,0,20000);
    TH2F *dopGenVQen = new TH2F("dopGenVQen","Doppler Corrected Gamma Energy vs Excitation Energy",10000,0,10000,20000,0,20000);
    TH2F *Crystals = new TH2F("Crystals","Gretina Crystals",360,0,360,180,0,180);
    TH2F *AddGenVQen = new TH2F("AddGenVQen","Add-back-ed Gamma Energy vs Excitation Energy",1000,0,10000,20000,0,20000);
    TH2F *dopAddGenVQen = new TH2F("dopAddGenVQen","Doppler Corrected Add-Back-ed Gamma Energy vs Excitation Energy",10000,0,10000,20000,0,20000);

		
    //ENERGY LOSS DEFINITIONS (DEADLAYER AND TARGET)	
	EnergyLoss* DeadLayer = new EnergyLoss("ProtonInSi.dat");
	EnergyLoss* Target = new EnergyLoss("ProtonInC2D4.dat");
	//DeadLayer->UseGL1024();

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
            std::cout << "Analyzing Run Number: " << runNumberStr << std::endl;
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
		

		//SX3 MULTIPLICITY LOOP STARTS HERE
		for(Int_t j=0; j<SX3Mul; j++){

            //Removing bad and emptry strips
			if (SX3Det[j]==5 || SX3Det[j]==11)continue;
            if (SX3Det[j]==2 && SX3Strip[j]==3)continue;
            if (SX3Det[j]==3 && SX3Strip[j]==3)continue;
            if (SX3Det[j]==6 && SX3Strip[j]==1)continue;
            if (SX3Det[j]==6 && SX3Strip[j]==2 && SX3Sector[j]==0)continue;
            if (SX3Det[j]==6 && SX3Strip[j]==3 && SX3Sector[j]==0)continue;
            if (SX3Det[j] > 11)continue; 

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
            if(SectorEnergy < 1e-6) continue;

			//IC Tracking Stuffs for reaction angle
			// Proton (x,y,z) for all SX3 detectors (proton hit position)
			Float_t radius = 98.;// Radial distance to all detectors from center of the barrel = 98mm
            Float_t SX3_Phi = SX3Det[j]*30. + ((3-SX3Strip[j]-1)*atan2(5.,radius)*180./M_PI) + ((3-SX3Strip[j]-2)*atan2(5.,radius)*180./M_PI);
            Float_t ProtonX = radius * sin (SX3_Phi*M_PI/180.); 
			Float_t ProtonY = radius * cos(SX3_Phi*M_PI/180.);
			Float_t ProtonZ = -1*(Position + 3.3);
			
			//Target (x,y,z) from Ionization Chamber using Scaling factors
            Float_t ScaleX = 0.3; 
			Float_t ScaleY = 0.36; 
			Float_t ShiftX = 0;
	        Float_t ShiftY = 0;

            Float_t TargetX = (icPositionX*ScaleX) + ShiftX;
            Float_t TargetY = (icPositionY*ScaleY) + ShiftY;
            Float_t TargetZ = TargetY*(tan(27.*M_PI/180.));

            //Defining two vectors to find angle
			TVector3 TargetToOriginVector(0,0,ProtonZ-TargetZ); //This vector Points from Target spot to Centre of QQQ5 Det
			TVector3 TargetToProtonVector(ProtonX-TargetX,ProtonY-TargetY,ProtonZ-TargetZ); // This vector points from Target spot to Proton spot in QQQ5.
			Float_t SX3Angle = 180. - ((TargetToOriginVector.Angle(TargetToProtonVector))*180./M_PI);
            
            //Effective DeadLayer and Target Thickness
            Float_t DeadLayer_thickness;
            targ_thick((180.-SX3Angle)*M_PI/180., SX3_Phi*M_PI/180., 0.0001, 0.0001, 90*M_PI/180., &DeadLayer_thickness);//Dead layer thickness = 0.0001mm and Depth = thickness
            Float_t Target_thickness;
            targ_thick((180.-SX3Angle)*M_PI/180., SX3_Phi*M_PI/180., 0.00292097835, 0.00292097835, 333*M_PI/180., &Target_thickness); //Target thickness = 0.005841956mm and depth = 0.5*target thickness
            
			//Energy Loss Calculations (DeadLayer and Target)
			Float_t SectorEnergyDeadLayer = DeadLayer->AddBack(SectorEnergy/1000., DeadLayer_thickness); 
			Double_t SectorEnergyTarget = 1000.*( Target->AddBack(SectorEnergyDeadLayer, 0.00292097835/fabs(cos((153.-SX3Angle)*M_PI/180.))) );  

			//Excitation Energy Calculation
			Double_t SectorExcitationEnergy = (Qgs) -  (1./My) * ( (My+Mb)*(SectorEnergyTarget/1000.) - (My-Ma)*Ta - (2.*sqrt(Ma*Mb*Ta*(SectorEnergyTarget/1000.))* cos( SX3Angle * (M_PI/180.) ) ) ) ;

            //Finding the beta for doppler correction
            Float_t QVal = Qgs - SectorExcitationEnergy;
	        Float_t Ty = QVal + Ta - (SectorEnergyTarget/1000.); //Recoil Energy (Q=Ty+Tb-Tx-Ta......From Krane pg.381)
	        Float_t beta = (1/(Ty+(My*931.5))) * (sqrt((2*Ty*My*931.5)+(Ty*Ty)));

            //Finding the phi angle using IC
            TVector3 Vec1(0,0,390);
            TVector3 Vec2(icPositionX-TargetX, icPositionY-TargetY,390);
            Float_t IC_phi = ( Vec1.Angle(Vec2) )*180./M_PI;

            //Finding the phi angle using kinematics formula
            Float_t Kin_phi = asin(sqrt((Mb*SectorEnergyTarget/1000.)/(My*Ty)) * sin((SX3Angle)*M_PI/180.) ) * 180./M_PI;

            //Finding x,y for recoil using the Kinematics phi
            Float_t Recoil_theta = Kin_phi;
            Float_t Recoil_phi = SX3_Phi + 180.;
            Float_t Recoil_R = 390./(cos(Recoil_theta*M_PI/180.));
            Float_t Recoil_X = Recoil_R * sin(Recoil_theta * M_PI/180.) * sin(Recoil_phi * M_PI/180.);
            Float_t Recoil_Y = Recoil_R * sin(Recoil_theta * M_PI/180.) * cos(Recoil_phi * M_PI/180.);
            Float_t Recoil_Z = 390.;


            //PROCESSING XTAL SINGLES
            //Variables Declaration
            Float_t X_xtals[xtalsMul]={0}, Y_xtals[xtalsMul]={0}, Z_xtals[xtalsMul]={0};
            Float_t Xtals_Theta[xtalsMul]={0}, Xtals_Phi[xtalsMul]={0};

            Float_t xtals_addback[(xtalsMul-1)*(xtalsMul-1)]={0};
            Float_t Xtals_addback_eff[(xtalsMul-1)*(xtalsMul-1)]={0};
            Float_t X_addback[(xtalsMul-1)*(xtalsMul-1)]={0}, Y_addback[(xtalsMul-1)*(xtalsMul-1)]={0}, Z_addback[(xtalsMul-1)*(xtalsMul-1)]={0};
            Int_t AddbackID[xtalsMul] = {0};

            //Loop over xtals Multiplicity
            for (Int_t k=0; k<xtalsMul; k++){

                //Offset correction and Transformation to make x,y,z consistent in ORRUBA and GRETINA
                X_xtals[k] = (-1)*xtals_ylab[k];   //Crystals_ylab is "-x" in my ORRUBA geometry
                Y_xtals[k] = (-1)*xtals_xlab[k];   //Crystals_xlab() is "-y" in my ORRUBA geometry
                Z_xtals[k] = xtals_zlab[k];        //Crystals_zlab() is "z" in my ORRUBA geometry

                //GRETINA OFFSET
                if(X_xtals[k]<0.){
                    X_xtals[k] = X_xtals[k] - GretinaOffsetX1; // Beam Left Shift
                } else{
                    X_xtals[k] = X_xtals[k] + GretinaOffsetX2; //Beam Right Shift
                }
                Y_xtals[k] = Y_xtals[k] - GretinaOffsetY; //ORRUBA  higher than GRETINA
                Z_xtals[k] = Z_xtals[k] - GretinaOffsetZ; //ORRUBA downstream of GRETINA

	            //Finding the gammaAngle for doppler correction
	            //If we assume pin-point beam at origin
	            TVector3 V2((X_xtals[k]),(Y_xtals[k]),(Z_xtals[k]));
	            TVector3 V3(0,0,390);
	            Float_t gAngle = (V3.Angle(V2))*180./M_PI;
                                
	            //Doppler Correcting here
                Float_t dGamma = dop(gAngle,beta,xtals_cc[k]);

                //Finding polar coordinates
                 Float_t tempf = sqrt( (pow(X_xtals[k],2.))+(pow(Y_xtals[k],2.))+(pow(Z_xtals[k],2.)) );
                Float_t tempf2 = sqrt( (pow(X_xtals[k],2.))+(pow(Y_xtals[k],2.)) );
                if(Z_xtals[k]/tempf >= -1. && Z_xtals[k]/tempf >= -1.){
                    Xtals_Theta[k] = acos(Z_xtals[k]/tempf);
                }else{
                    std::cout << "cos(gam_theta) = " << Z_xtals[k]/tempf << " was out of range" << std::endl;
                    std::cout << "X_xtals = " << X_xtals[k] << " Y_xtals = " << Y_xtals[k] << "Z_xtals = " << Z_xtals[k] << std::endl; 
                    Xtals_Theta[k] = 0;
                }
                if(Y_xtals[k]/tempf2 >=-1. && Y_xtals[k]/tempf2 <=1.){ //Mithch's X is Y in my Coordinate sys
                    Xtals_Phi[k] = acos(Y_xtals[k]/tempf2);
                    if(X_xtals[k]<0)Xtals_Phi[k] = 2.0*M_PI - Xtals_Phi[k];
                }else{
                    std::cout << "cos(gam_phi) = " << Xtals_Phi[k] << " was out of range" << std::endl;
                    std::cout << "X_xtals = " << X_xtals[k] << " Y_xtals = " << Y_xtals[k] << "Z_xtals = " << Z_xtals[k] << std::endl; 
                    Xtals_Phi[k] = 0;
                }

                //Filling Gamma histograms
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950.) GenVQen->Fill(xtals_cc[k],(SectorExcitationEnergy*1000.));
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950.) dopGenVQen->Fill(dGamma,(SectorExcitationEnergy*1000.));
				if (SX3Upstream[j] && tdcIC>600. && tdcIC<950.) Crystals->Fill(Xtals_Phi[k]*180./M_PI,Xtals_Theta[k]*180./M_PI);
            }//End of XTAL Singles    

			
			//PROCESSING XTAL ADDBACK
            //"Cluster" addback
            Int_t gmult=0;
            for (Int_t k=0; k<xtalsMul; k++){
                for (Int_t l=0; l<xtalsMul; l++){
                    if(k!=l && AddbackID[k]==0 && AddbackID[l]==0){
                        Float_t tempf = sin(Xtals_Theta[k])*sin(Xtals_Theta[l])*cos(Xtals_Phi[k]-Xtals_Phi[l]) + cos(Xtals_Theta[k])*cos(Xtals_Theta[l]);
                        if (tempf >= cos(ClusterAngle*M_PI/180.) && tempf <= cos(0.0)){

                            //switch the addback id
                            AddbackID[k]=1;
                            AddbackID[l]=1;

                            //Adding the Gamma energies
                            xtals_addback[gmult] = xtals_cc[k] + xtals_cc[l];                          
                            //Choosing (x,y,z) for doppler correction 
                            if (xtals_cc[k] > xtals_cc[l]){
                                X_addback[gmult] = X_xtals[k];
                                Y_addback[gmult] = Y_xtals[k];
                                Z_addback[gmult] = Z_xtals[k];
                            }else{
                                X_addback[gmult] = X_xtals[l];
                                Y_addback[gmult] = Y_xtals[l];
                                Z_addback[gmult] = Z_xtals[l];
                            }
                            gmult++;
                        }
                    }
                }//End of loop over l  
                //Keeping the gammas that are not added back
                if(AddbackID[k]==0){
                    xtals_addback[gmult] = xtals_cc[k];
                    X_addback[gmult] = X_xtals[k];
                    Y_addback[gmult] = Y_xtals[k];
                    Z_addback[gmult] = Z_xtals[k];
                    gmult++;
                }
            }//END OF ADDBACK   

			//DOPPLER CORRECTION ON ADD-BACK-ED GAMMAS
            Float_t addback_gAngle_Origin[gmult];
            for(Int_t k=0; k<gmult; k++){
                //If we assume pin-point beam at origin
	            TVector3 V2((X_addback[k]),(Y_addback[k]),(Z_addback[k]));
	            TVector3 V3(0,0,390);
	            Float_t addback_gAngle = (V3.Angle(V2))*180./M_PI;
                
	            //Doppler Correcting here
                Float_t addback_dGamma = dop(addback_gAngle,beta,xtals_addback[k]);
               
                //Fill your histogram after addback
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950.) AddGenVQen->Fill(xtals_addback[k],(SectorExcitationEnergy*1000.));
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950.) dopAddGenVQen->Fill(addback_dGamma,(SectorExcitationEnergy*1000.));    
            }//END 0F ADDBACK 

            //Fill your SX3 histograms here
			if(SX3Upstream[j])SX3EvA->Fill(SX3Angle,SectorEnergyTarget);
			if(SX3Upstream[j] && tdcIC>600. && tdcIC<950.)SX3EvA_TDC->Fill(SX3Angle,SectorEnergyTarget);
			if(SX3Upstream[j])SX3Ex->Fill(SectorExcitationEnergy);
			if(SX3Upstream[j] && tdcIC>600. && tdcIC<950.)SX3Ex_TDC->Fill(SectorExcitationEnergy);
        }//End of Loop over SX3 Multiplicity    
    }// End of event by event analysis
  
    
	//WRITING HISTOGRAM
    outputFile->cd();

	//IC Stuffs
    hICdEE->Write();
    hICdEE_cut->Write();
    
    //SX3 Stuffs
    SX3EvA->Write();
    SX3EvA_TDC->Write();
    SX3Ex->Write();
    SX3Ex_TDC->Write();

	//Gamma Stuffs
    GenVQen->Write();
    dopGenVQen->Write();
    Crystals->Write();
    AddGenVQen->Write();
    dopAddGenVQen->Write();

    outputFile->Close();
} //END OF THE ANALYSIS LOOP
