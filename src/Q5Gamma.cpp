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

//MAIN ANALYSIS LOOP
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
    TFile* outputFile = new TFile(OutputPath + OutputFilePrefix + "_QQQ5.root", "recreate");

    //OPEN ALL INPUT FILES
    //Open QQQ5 Energy Calibration File
	std::ifstream QQQ5EnCalFile;
	QQQ5EnCalFile.open("Calibration/UpQQQ5EnCal.dat");
	Float_t QQQ5EnCalSlope[512] = {0};
	Float_t QQQ5EnCalIntercept[512] = {0};
    for (Int_t i=0; i<512; i++){
        QQQ5EnCalFile >> QQQ5EnCalSlope[i] >> QQQ5EnCalIntercept[i];
        //std::cout << i << '\t' << QQQ5EnCalSlope[i] << '\t' << QQQ5EnCalIntercept[i] << std::endl;
    }
    QQQ5EnCalFile.close();

	//Open QQQ5 Angle file for Exciation Energy vs Angle bins
	std::ifstream afile;
	afile.open("Geometry/QQQ5ExcitationAngle.dat");
	Double_t Q5ExAngle[33] = {0};
	for (Int_t i=0; i<33; i++){
		afile >> Q5ExAngle[i];
		//std::cout << Q5ExAngle[i] << std::endl;
	}
	afile.close();

	//Open QQQ5 Geometry File for Ring distance for IC Tracking
	std::ifstream geofile;
	geofile.open("Geometry/QQQ5Geometry.dat");
	Double_t QQQ5RingDistance[32] = {0};
	Double_t QQQ5RingAngle[32]={0};
	for (Int_t i=0; i<32; i++){
		geofile >> QQQ5RingDistance[i] >> QQQ5RingAngle[i];
		//std::cout << QQQ5RingDistance[i] << std::endl;
	}
	geofile.close();

	//DEFINE HISTOGRAMS HERE
	//IC stuffs
	TH2F* hICdEE = new TH2F("icdEE", "Ionization Chamber dE vs E; E; dE", 500, 0, 4000, 500, 0, 4000);
    TH2F* hICdEE_cut = new TH2F("icdEE_cut", "Ionization Chamber dE vs E (with IC Cut); E; dE", 500, 0, 4000, 500, 0, 4000);	

	//QQQ5 Stuffs
	TH2F *Q5EvA = new TH2F("Q5EvA","QQQ5 Energy vs Angle",180,0,180,2000,0,20000);
    TH2F *Q5EvA_TDC = new TH2F("Q5EvA_TDC","QQQ5 Energy vs Angle with TDC",180,0,180,2000,0,20000);
	TH1F *Q5Ex = new TH1F("Q5Ex","QQQ5 Excitation Energy",250,-2,8);
	TH1F *Q5Ex_TDC = new TH1F("Q5Ex_TDC","QQQ5 Excitation Energy with TDC",250,-2,8);

	
    //GammaRay Stuffs
	TH2F *GenVQen = new TH2F("GenVQen","Gamma Energy vs Excitation Energy",1000,0,10000,20000,0,20000);
    TH2F *dopGenVQen = new TH2F("dopGenVQen","Doppler Corrected Gamma Energy vs Excitation Energy",10000,0,10000,20000,0,20000);
    TH2F *Crystals = new TH2F("Crystals","Gretina Crystals",360,0,360,180,0,180);
    TH2F *AddGenVQen = new TH2F("AddGenVQen","Add-back-ed Gamma Energy vs Excitation Energy",1000,0,10000,20000,0,20000);
    TH2F *dopAddGenVQen = new TH2F("dopAddGenVQen","Doppler Corrected Add-Back-ed Gamma Energy vs Excitation Energy",10000,0,10000,20000,0,20000);

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
        //Everything after this has the ic cut on dE vs E
	    hICdEE_cut->Fill(icE, icdE);

        //QQQ5 MULTIPLICITY LOOP STARTS HERE
	    for(Int_t j=0; j<QQQ5Mul; j++){
	        if(QQQ5RingEnergy[j] < 1e-6) continue;

	        //Proton x,y,z for QQQ5 detectors
	        Float_t Q5_Phi = 11.25 + (QQQ5Sector[j]*22.5) + (QQQ5Det[j]*90.);
	        Float_t ProtonX = QQQ5RingDistance[QQQ5Ring[j]]* sin(Q5_Phi*(M_PI/180.));
			Float_t ProtonY = QQQ5RingDistance[QQQ5Ring[j]]* cos(Q5_Phi*(M_PI/180.));
			Float_t ProtonZ = -80.; //Distance from Centre of target to QQQ5 Det

            //SCALING AND SHIFT FACTORS
            Float_t ScaleX[4] = {0,0,0,0}; 
			Float_t ScaleY[4] = {0,0,0,0};
			Float_t ShiftX[4] = {0,0,0,0};
			Float_t ShiftY[4] = {0,0,0,0};

            Float_t TargetX = (icPositionX*ScaleX[QQQ5Det[j]]) + ShiftX[QQQ5Det[j]];
			Float_t TargetY = (icPositionY*ScaleY[QQQ5Det[j]]) + ShiftY[QQQ5Det[j]];
			Float_t TargetZ = TargetY*tan(27.*M_PI/180.);

	        //Defining two vectors to find angle
			TVector3 TargetToOriginVector(0,0,ProtonZ-TargetZ); //This vector Points from Target spot to Centre of QQQ5 Det
			TVector3 TargetToProtonVector(ProtonX-TargetX,ProtonY-TargetY,ProtonZ-TargetZ); // This vector points from Target spot to Proton spot in QQQ5.
			Float_t Q5Angle = 180. - ((TargetToOriginVector.Angle(TargetToProtonVector))*180./M_PI);

            //Effective DeadLayer and Target Thickness
            Float_t DeadLayer_thickness;
            targ_thick(Q5Angle*M_PI/180., Q5_Phi*M_PI/180., 0.0001, 0.0001, 0*M_PI/180., &DeadLayer_thickness);//Dead layer thickness = 0.0001mm and Depth = thickness
            Float_t Target_thickness;
            targ_thick((180.-Q5Angle)*M_PI/180., Q5_Phi*M_PI/180., 0.00292097835, 0.00292097835, 27*M_PI/180., &Target_thickness); //Target thickness = 0.005841956mm and depth = 0.5*target thickness
           
            //Energy Calibration
            Double_t QQQ5EnergyRing = QQQ5RingADC[j]*QQQ5EnCalSlope[QQQ5Ring[j]+(QQQ5Sector[j]*32)+(QQQ5Det[j]*128)] + QQQ5EnCalIntercept[QQQ5Ring[j]+(QQQ5Sector[j]*32)+(QQQ5Det[j]*128)];
			if(QQQ5EnergyRing < 1e-6) continue;

	        //Energy Loss Calculation 
            Double_t QQQ5RingEnergyDeadLayer = DeadLayer->AddBack(QQQ5EnergyRing/1000., DeadLayer_thickness);
	        Double_t QQQ5RingEnergyTarget = Target->AddBack(QQQ5RingEnergyDeadLayer, Target_thickness); 
            
	        //Excitation Energy and Q-Value Calculations
	        Double_t RingExcitationEnergy = (Qgs) -  (1./My) * ( (My+Mb)*QQQ5RingEnergyTarget - (My-Ma)*Ta - (2.*sqrt(Ma*Mb*Ta*QQQ5RingEnergyTarget)* cos( Q5Angle * (M_PI/180.) ) ) ) ;
	        Double_t RingQValue = (1./My) * ( (My+Mb)*QQQ5RingEnergyTarget - (My-Ma)*Ta - (2.*sqrt(Ma*Mb*Ta*QQQ5RingEnergyTarget)* cos( Q5Angle * (M_PI/180.) ) ) ) ;
			
	        //Finding the beta for doppler correction
	        Float_t Ty = RingQValue + Ta - QQQ5RingEnergyTarget; //Recoil Energy
            Float_t beta =  (sqrt((2*Ty*My*931.494013)+(Ty*Ty))) / (Ty+(My*931.494013));

            //Finding the phi angle using IC
            TVector3 Vec1(0,0,390);
            TVector3 Vec2(icPositionX-TargetX, icPositionY-TargetY,390);
            Float_t IC_phi = ( Vec1.Angle(Vec2) )*180./M_PI;

            //Finding the phi angle using kinematics formula
            Float_t Kin_phi = asin(sqrt((Mb*QQQ5RingEnergyTarget)/(My*Ty)) * sin((180.-Q5Angle)*M_PI/180.) ) * 180./M_PI; //This is actually theta angle for recoil

            //Finding x,y for recoil using the Kinematics phi
            Float_t Recoil_theta = Kin_phi;
            Float_t Recoil_phi = Q5_Phi + 180.;
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

                //Fill your Gamma Histograms Here
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) GenVQen->Fill(xtals_cc[k],(RingExcitationEnergy*1000.));
	            if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) dopGenVQen->Fill(dGamma,(RingExcitationEnergy*1000.));  
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) Crystals->Fill(Xtals_Phi[k]*180./M_PI,Xtals_Theta[k]*180./M_PI);  
 
            }//End of XTAL Singles    

            //XTAL ADDBACK STARTS HERE
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

            //DOPPLER CORRECTION ON ADD-BACKE-ED GAMMAS
            Float_t addback_gAngle_Origin[gmult];
            for(Int_t k=0; k<gmult; k++){
                //If we assume pin-point beam at origin
	            TVector3 V2((X_addback[k]),(Y_addback[k]),(Z_addback[k]));
	            TVector3 V3(0,0,390);
	            Float_t addback_gAngle = (V3.Angle(V2))*180./M_PI;
                
	            //Doppler Correcting here
                Float_t addback_dGamma = dop(addback_gAngle,beta,xtals_addback[k]);

                //Fill your histogram after addback
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) AddGenVQen->Fill(xtals_addback[k],(RingExcitationEnergy*1000.));     
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) dopAddGenVQen->Fill(addback_dGamma,(RingExcitationEnergy*1000.));     
            }//END OF ADDBACK

	        //FILL YOUR HISTOGRAMS HERE
	        // QQQ5 Stuffs
            if(QQQ5Upstream[j]) Q5Ex->Fill(RingExcitationEnergy);
	        if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) Q5Ex_TDC->Fill(RingExcitationEnergy);
	        if(QQQ5Upstream[j])Q5EvA->Fill(Q5Angle,(1000.*QQQ5RingEnergyTarget));
            if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.)Q5EvA_TDC->Fill(Q5Angle,(1000.*QQQ5RingEnergyTarget));
        }//End of QQQ5 Multiplicity
    }// End of event by event analysis

    //WRITING HISTOGRAM
    outputFile->cd();

    //WRITING HISTOGRAMS TO OUTPUT FILE HERE
    //IC stuffs
    hICdEE->Write();
    hICdEE_cut->Write();
    
    //QQQ5 Stuffs
    Q5EvA->Write();
    Q5EvA_TDC->Write();
    Q5Ex->Write();
    Q5Ex_TDC->Write();

    //Gamma Stuffs
    GenVQen->Write();
    dopGenVQen->Write();
    Crystals->Write();
    AddGenVQen->Write();
    dopAddGenVQen->Write();
    
    outputFile->Close();

} //END OF THE ANALYSIS LOOP

