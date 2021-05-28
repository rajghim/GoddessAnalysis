#define Analysis_cxx
#include "Analysis.h"
#include "EnergyLoss.h"
#include <TRandom3.h>
#include <string>
#include "EffectiveThickness.h"

//----------------------------------------------------------------Doppler Correction Function Definition---------------------------------------------------------------------//
Float_t dop(Float_t GamAngle, Float_t RecoilBeta, Float_t GamEnergy){
    Float_t RecoilGamma = 1 / (sqrt (1-(RecoilBeta*RecoilBeta) ) );
    Float_t DopplerCorrectedEnergy = GamEnergy * ( 1-RecoilBeta*(cos(GamAngle*M_PI/180.)) ) * RecoilGamma;
    return DopplerCorrectedEnergy;
}

// -----------------------------------------------------------Main Analysis Loop-------------------------------------------------------------------------------------------//

void Analysis::Loop() {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    TRandom3 *r = new TRandom3();

    TString CutPath = "/mnt/e/Analysis/Analysis Code/Cuts/ICCuts/";
    TString CutPrefix = "Pcut_p30dp";

    int prevRunNumber = -1;

//--------------------------------------------------------------Open all input files-----------------------------------------------------------------------------------------//

    //Open the pedestals File (Keep Left Pedestals and Right Pedestals and Back Pedestals in the same file) 
	std::ifstream file;
	file.open("Calibration/SX3pedestals.dat");
	Double_t LeftPed[192] = {0};
	Double_t RightPed[192] = {0};
	Double_t BackPed[192] = {0};
	for (Int_t i = 0; i<192;i++){
		file >> LeftPed[i] >> RightPed[i] >> BackPed[i];	
		//std::cout << LeftPed[1] << '\t' << RightPed[5] << '\t' << BackPed[50] << std::endl;
	}
 
	//Open the gains file
	std::ifstream gainfile;
	gainfile.open("Calibration/SX3gains.dat");
	Double_t Gains[48] = {0}; 
	for (Int_t i=0; i<48; i++){
			gainfile >>Gains[i];
			//std::cout << Gains[2] << std::endl;
	}



	//Open the Energy Calibration file for Upstream SX3 (Keep both Front and Back Encal Parameters in same file)
	std::ifstream SX3EnCalfile;
	SX3EnCalfile.open("Calibration/SX3EnCal.dat");
	Double_t SX3EnCalSlope[192] = {0};
	Double_t SX3EnCalIntercept[192] = {0};
	Double_t SX3SectorEnCalSlope[192] = {0};
	Double_t SX3SectorEnCalIntercept[192] = {0}; 
	for (Int_t i=0; i<192; i++){
			SX3EnCalfile >> SX3EnCalSlope[i] >> SX3EnCalIntercept[i] >> SX3SectorEnCalSlope[i] >> SX3SectorEnCalIntercept[i];
	}

	//Open the Position Calibration file for Upstream SX3
	std::ifstream SX3PosCalfile;
	SX3PosCalfile.open("Calibration/SX3PosCal.dat");
	Double_t xone[48] = {0};
	Double_t xtwo[48] = {0}; 
	for (Int_t i=0; i<48; i++){
			SX3PosCalfile >> xone[i] >> xtwo[i];
	}


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

	//Open QQQ5 Geometry File for Ring distance for IC Tracking
	std::ifstream roughgeofile;
	roughgeofile.open("Geometry/roughQQQ5Geometry.dat");
	Double_t roughQQQ5InRingDistance[32] = {0};
	Double_t roughQQQ5OutRingDistance[32] = {0};
	for (Int_t i=0; i<32; i++){
		roughgeofile >> roughQQQ5InRingDistance[i] >> roughQQQ5OutRingDistance[i];
		//std::cout << roughQQQ5InRingDistance[i] << '\t' << roughQQQ5OutRingDistance[i] << std::endl;
	}
	geofile.close();

	//QuadCutoff file
	std::ifstream quadcutofffile;
	quadcutofffile.open("Geometry/QuadCutoff.dat");
	Double_t quadcut[11] = {0};
	for (Int_t i=0; i<11; i++){
		quadcutofffile >> quadcut[i];
		//std::cout << quadcut[i] << std::endl;
	}
	quadcutofffile.close();

    //Gains and Offsets to fix the shifts in excitation energy plots of QQQ5 and SX3 
    Float_t Q5ExEnGain[4] = {1.002322,1.005295,1.004755,0.987985};
	Float_t Q5ExEnOffset[4] = {-0.06254,-0.2584,-0.24833,-0.20053};
    Float_t SX3ExEnGain[12] = {0.934352145,1.018363202,0.966706885,0.962783882,0.900456366,1,0.974268997,1.00211885,1.003156092,0.967291741,1.018601585,1};
	Float_t SX3ExEnOffset[12] = {-0.246132018,-0.391123183,-0.270966644,-0.343275319,-0.111066209,0,-0.31608631,-0.401416842,-0.413230614,-0.386436621,-0.295758643,0};


//----------------------------------------------------------------------Histogram Definition--------------------------------------------------------------------------------------//

	//Define Histograms here
	//IC stuffs
	TH2F* hICdEE = new TH2F("icdEE", "Ionization Chamber dE vs E; E; dE", 500, 0, 4000, 500, 0, 4000);
    TH2F* hICdEE_cut = new TH2F("icdEE_cut", "Ionization Chamber dE vs E (with IC Cut); E; dE", 500, 0, 4000, 500, 0, 4000);
    TH2F* hICdEE_kincut = new TH2F("icdEE_kincut", "Ionization Chamber dE vs E (with IC Cut and (d,p) kinematics cut); E; dE", 500, 0, 4000, 500, 0, 4000);
	
	//QQQ5 Stuffs
	TH2F *Q5_RvS_TDC = new TH2F("Q5_RvS_TDC","Ring Vs Sector Energy for QQQ5 with TDC",2000,0,20,2000,0,20);
    TH2F *Q5_RvS_TDC_cut = new TH2F("Q5_RvS_TDC_cut","Ring Vs Sector Energy for QQQ5 with TDC with the cut",2000,0,20,2000,0,20);
	TH2F *Q5EvA_TDC = new TH2F("Q5EvA_TDC","QQQ5 Energy vs Angle with TDC",1000,0,180,2000,0,20000);
    TH2F *Q5EvA_TDC_kinCut = new TH2F("Q5EvA_newTDC_kinCut","QQQ5 Energy vs Angle with new TDC and cut on (d,p) events",1000,0,180,2000,0,20000);
	TH1F *Q5Ex_TDC = new TH1F("Q5Ex_TDC","QQQ5 Excitation Energy with TDC",250,-2,8);
	
    //SX3 Stuffs
    TH1F *SX3Ex_TDC = new TH1F("SX3Ex_TDC","SX3 Excitation Energy with TDC",250,-2,8);
    TH2F *SX3_FvB_TDC = new TH2F("SX3_FvB_TDC","Sector vs Strip Energy for SX3 with TDC",2000,0,20000,2000,0,20000);
    TH2F *SX3_FvB_TDC_cut = new TH2F("SX3_FvB_TDC_cut","Sector vs Strip Energy for SX3 with TDC and cut",2000,0,20000,2000,0,20000);
	
    //GammaRay Stuffs
	TH2F *GenVQen = new TH2F("GenVQen","Gamma Energy vs Excitation Energy",1000,0,10000,20000,0,20000);
    TH2F *GenVQen_Q5 = new TH2F("GenVQen_Q5","Gamma Energy vs QQQ5 Excitation Energy",1000,0,10000,20000,0,20000);
	TH2F *GenVQen_SX3 = new TH2F("GenVQen_SX3","Gamma Energy vs SX3 Excitation Energy",1000,0,10000,20000,0,20000);
    TH2F *dopGenVQen_Timecut = new TH2F("dopGenVQen_Timecut","Doppler Corrected Gamma Energy vs Excitation Energy with timestamp cut",1000,0,10000,20000,0,20000);
    TH2F *dopGenVQen_Timecut_Q5 = new TH2F("dopGenVQen_Timecut_Q5","Doppler Corrected Gamma Energy vs QQQQ5 Excitation Energy with timestamp cut",1000,0,10000,20000,0,20000);
    TH2F *dopGenVQen_Timecut_SX3 = new TH2F("dopGenVQen_Timecut_SX3","Doppler Corrected Gamma Energy vs SX3 Excitation Energy with timestamp cut",1000,0,10000,20000,0,20000);
    //Addback Histograms
    TH2F *dopAddGenVQen_Timecut = new TH2F("dopAddGenVQen_Timecut","Added Back & Doppler Corrected Gamma Energy vs Excitation Energy with timestamp cut",1000,0,10000,20000,0,20000);
    TH2F *dopAddGenVQen_Timecut_Q5 = new TH2F("dopAddGenVQen_Timecut_Q5","Added Back & Doppler Corrected Gamma Energy vs QQQQ5 Excitation Energy with timestamp cut",1000,0,10000,20000,0,20000);
    TH2F *dopAddGenVQen_Timecut_SX3 = new TH2F("dopAddGenVQen_Timecut_SX3","Added Back & Doppler Corrected Gamma Energy vs SX3 Excitation Energy with timestamp cut",1000,0,10000,20000,0,20000);
    
    //Particle Histograms With Cuts in Gammas
    TH1F *Q5Ex_TDC_1266 = new TH1F("Q5Ex_TDC_1266","QQQ5 Excitation Energy with TDC and 1266 keV.13 Gamma Ray (1266.13->0)",250,-2,8);
    TH1F *SX3Ex_TDC_1266 = new TH1F("SX3Ex_TDC_1266","SX3 Excitation Energy with TDC and 1266.13 keV Gamma Ray (1266.13->0)",250,-2,8);
    TH1F *Q5Ex_TDC_3134 = new TH1F("Q5Ex_TDC_3134","QQQ5 Excitation Energy with TDC and 3134.3 keV Gamma Ray (3134.3->0)",250,-2,8);
    TH1F *SX3Ex_TDC_3134 = new TH1F("SX3Ex_TDC_3134","SX3 Excitation Energy with TDC and 3134.3 keV Gamma Ray (3134.3->0)",250,-2,8);  
    TH1F *Q5Ex_TDC_2028 = new TH1F("Q5Ex_TDC_2028","QQQ5 Excitation Energy with TDC and 2028.8 keV Gamma Ray (3295->1266.13)",250,-2,8);
    TH1F *SX3Ex_TDC_2028 = new TH1F("SX3Ex_TDC_2028","SX3 Excitation Energy with TDC and 2028.8 keV Gamma Ray (3295->1266.13)",250,-2,8);
    TH1F *Q5Ex_TDC_2148 = new TH1F("Q5Ex_TDC_2148","QQQ5 Excitation Energy with TDC and 2148.4 keV Gamma Ray (3414.6->1266.13)",250,-2,8);
    TH1F *SX3Ex_TDC_2148 = new TH1F("SX3Ex_TDC_2148","SX3 Excitation Energy with TDC and 2148.4 keV Gamma Ray (3414.6->1266.13)",250,-2,8);
    TH1F *Q5Ex_TDC_3506 = new TH1F("Q5Ex_TDC_3506","QQQ5 Excitation Energy with TDC and 3506 keV Gamma Ray (3506->0)",250,-2,8);
    TH1F *SX3Ex_TDC_3506 = new TH1F("SX3Ex_TDC_3506","SX3 Excitation Energy with TDC and 3506 keV Gamma Ray (3506->0)",250,-2,8);

    

//--------------------------------------------------------------------Reaction Paramteters, Output files and Cut files------------------------------------------------------//

	//Create Output File   
	TFile* outputFile = new TFile("/mnt/e/Analysis/Analysis Code/Output/Goddess30P.root", "recreate");
	
	//Reaction Parameters
	Float_t Ma = 29.97831; //Mass of beam (30P)
	Float_t Mx = 2.014102; //Mass of target (deuterium)
	Float_t Mb = 1.007276466; //Mass of Ejectile (proton) 
	Float_t My = 30.97376199; //Mass of Recoil (31P)
	Float_t Qgs = 10.086439; //Q-value of Ground State
	Float_t Ta = 240.298; //Beam Energy(After the energy loss in Target)

    //Define Gretina Offsets Here
    Float_t GretinaOffset_X1 = 7.061; //Beam Left
    Float_t GretinaOffset_X2 = 12.37; //Beam Right
    Float_t GretinaOffset_Y = 0.13; //ORRUBA 0.13mm higher than GRETINA
    Float_t GretinaOffset_Z = 0.2325; //ORRUBA 0.2325mm downstream than GRETINA

    Float_t ClusterAngle = 15.;

	EnergyLoss* DeadLayer = new EnergyLoss("ProtonInSi.dat");
	EnergyLoss* Target = new EnergyLoss("ProtonInC2D4.dat");
	//DeadLayer->UseGL1024();

    //Front vs Back QQQ5 cut
    TCutG *QQQ5Cut_FvB;
    TFile *QQQ5CutFile;
    TString QQQ5CutFileName = "/mnt/e/Analysis/Analysis Code/Cuts/FvBcut_QQQ5.root";
    QQQ5CutFile = TFile::Open(QQQ5CutFileName);
    TString QQQ5CutName = "q5Cut";
    QQQ5Cut_FvB = static_cast<TCutG*>(QQQ5CutFile->Get(QQQ5CutName));
    QQQ5CutFile->Close();

    //Front vs Back SX3 cut
    TCutG *SX3Cut_FvB;
    TFile *SX3CutFile;
    TString SX3CutFileName = "/mnt/e/Analysis/Analysis Code/Cuts/FvBcut_SX3.root";
    SX3CutFile = TFile::Open(SX3CutFileName);
    TString SX3CutName = "sx3Cut";
    SX3Cut_FvB = static_cast<TCutG*>(SX3CutFile->Get(SX3CutName));
    SX3CutFile->Close();

    //kinematics cut for (d,p) events from kinematics plot in upstream channels.
    TCutG *kinCut;
    TFile *kinCutFile;
    TString kinCutFileName = "/mnt/e/Analysis/Analysis Code/Cuts/kinematicsCut.root";
    kinCutFile = TFile::Open(kinCutFileName);
    TString kinCutName = "kincut";
    kinCut = static_cast<TCutG*>(kinCutFile->Get(kinCutName));
    kinCutFile->Close();


//-----------------------------------------------------------------Event by Event Analysis Start-----------------------------------------------------------------------------//

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
            TString cutFileName = CutPath + CutPrefix + runNumberStr + ".root";
            std::cout << "Analyzing Run Number: " << runNumberStr << std::endl;
            cutFile = TFile::Open(cutFileName);
            TString cutName = Form("Pcut%d", runNumber);
            icdEECut = static_cast<TCutG*>(cutFile->Get(cutName));
            prevRunNumber = runNumber;
            cutFile->Close();
        }
        hICdEE->Fill(icE, icdE);

        //Applying the IC cut here
        if(!icdEECut->IsInside(icE, icdE)) continue; //This was initial IC Cut

        // Everything after this has the ic cut on dE vs E
	    hICdEE_cut->Fill(icE, icdE);

//----------------------------------------------------------QQQ5 Multiplicity Loop start---------------------------------------------------------------------------------------------//        
	    for(Int_t j=0; j<QQQ5Mul; j++){
	        if(QQQ5RingEnergy[j] < 1e-6) continue;

	        //IC Tracking
	        //Proton x,y,z for QQQ5 detectors
	        Float_t Q5_Phi = 11.25 + (QQQ5Sector[j]*22.5) + (QQQ5Det[j]*90.);
	        Float_t ProtonX = QQQ5RingDistance[QQQ5Ring[j]]* sin(Q5_Phi*(M_PI/180.));
			Float_t ProtonY = QQQ5RingDistance[QQQ5Ring[j]]* cos(Q5_Phi*(M_PI/180.));
			Float_t ProtonZ = -80.; //Distance from Centre of target to QQQ5 Det

            //Scaling and Shift factors
            //Float_t ScaleX[4] = {0.14,0.5,0.2,0.08}; 
			//Float_t ScaleY[4] = {0.08,0,0.1,0.5};
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
           
	        //Energy Loss Calculation 
            Double_t QQQ5RingEnergyDeadLayer = DeadLayer->AddBack(QQQ5RingEnergy[j]/1000., DeadLayer_thickness);
	        Double_t QQQ5RingEnergyTarget = Target->AddBack(QQQ5RingEnergyDeadLayer, Target_thickness); 
		
	        //Excitation Energy and Q-Value Calculations
	        Double_t RingExcitationEnergy = (Qgs) -  (1./My) * ( (My+Mb)*QQQ5RingEnergyTarget - (My-Ma)*Ta - (2.*sqrt(Ma*Mb*Ta*QQQ5RingEnergyTarget)* cos( Q5Angle * (M_PI/180.) ) ) ) ;
            RingExcitationEnergy = RingExcitationEnergy*Q5ExEnGain[QQQ5Det[j]] + Q5ExEnOffset[QQQ5Det[j]];

	        Double_t RingQValue = (1./My) * ( (My+Mb)*QQQ5RingEnergyTarget - (My-Ma)*Ta - (2.*sqrt(Ma*Mb*Ta*QQQ5RingEnergyTarget)* cos( Q5Angle * (M_PI/180.) ) ) ) ;
			
	        //Finding the beta and gamma for doppler correction
	        Float_t Ty = RingQValue + Ta - QQQ5RingEnergyTarget; //Recoil Energy (Q=Ty+Tb-Tx-Ta......From Krane pg.381)
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


//---------------------------------------------------------------------Crystal Multiplicity Loop Start-------------------------------------------------------------------------------//
           //Xtal Singles
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
                    X_xtals[k] = X_xtals[k] - GretinaOffset_X1; // Beam Left Shift
                } else{
                    X_xtals[k] = X_xtals[k] + GretinaOffset_X2; //Beam Right Shift
                }
                Y_xtals[k] = Y_xtals[k] - GretinaOffset_Y; //ORRUBA  higher than GRETINA
                Z_xtals[k] = Z_xtals[k] - GretinaOffset_Z; //ORRUBA downstream of GRETINA

	            //Finding the gammaAngle for doppler correction
	            //If we assume pin-point beam at origin
                TVector3 V1((X_xtals[k]-TargetX),(Y_xtals[k]-TargetY),(Z_xtals[k]-TargetZ)); 
	            TVector3 V2((X_xtals[k]),(Y_xtals[k]),(Z_xtals[k]));
	            TVector3 V3(0,0,390);
	            Float_t gAngle_Origin = (V3.Angle(V2))*180./M_PI;
	
	            //This assumes that recoil hits the (icPositionX,icPositionY,390) in IC OR (Recoil_X,Recoil_Y,390) from calculation (CHOOSE ONE)
	            TVector3 V4(icPositionX-TargetX, icPositionY-TargetY, 390-TargetZ); //This uses (icPositionX,icPositionY,390)
                //TVector3 V4(Recoil_X-TargetX, Recoil_Y-TargetY, 390-TargetZ); //This uses (Recoil_X,Recoil_Y,390)
	            Float_t gAngle_IC = (V4.Angle(V1))*180./M_PI;
                
	            //Doppler Correcting here
                Float_t dGamma_Origin = dop(gAngle_Origin,beta,xtals_cc[k]);
                Float_t dGamma_IC = dop(gAngle_IC,beta,xtals_cc[k]);


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
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) GenVQen_Q5->Fill(xtals_cc[k],(RingExcitationEnergy*1000.));

                //Using other cuts here
                //Front vs Back cut in QQQ5 Detector
                if(!QQQ5Cut_FvB->IsInside(QQQ5RingEnergy[j]/1000.,QQQ5SectorEnergy[j]/1000.)) continue;
                //(d,p) events cut in kinematics
                //if(!kinCut->IsInside(Q5Angle, (QQQ5RingEnergyTarget*1000))) continue;

                //Filling the histogram after using the Front vs Back cut in QQQ5 and the (d,p) kinematics cut
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && timeStamp>0 && GRETINATimeStamp>0 && foundGRETINA && (timeStamp-GRETINATimeStamp)>20. && (timeStamp-GRETINATimeStamp)<60.) dopGenVQen_Timecut->Fill(dGamma_Origin,(RingExcitationEnergy*1000.));   
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && timeStamp>0 && GRETINATimeStamp>0 && foundGRETINA && (timeStamp-GRETINATimeStamp)>20. && (timeStamp-GRETINATimeStamp)<60.) dopGenVQen_Timecut_Q5->Fill(dGamma_Origin,(RingExcitationEnergy*1000.));   
            }//End of XTAL Singles    

//---------------------------------------------------------- XTAL ADDBACK -------------------------------------------------------------------------------------//	       
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
            }//End of loop over k   

//--------------------------------------------------------------XTAL ADDBACK ENDS & DOPPLER CORRECTING "ADDBACK-ED" GAMMAS STARTS HERE--------------------------------------------------------------------------------------------------//
            Float_t addback_gAngle_Origin[gmult];
            for(Int_t k=0; k<gmult; k++){
                //If we assume pin-point beam at origin
                TVector3 V1((X_addback[k]-TargetX),(Y_addback[k]-TargetY),(Z_addback[k]-TargetZ)); 
	            TVector3 V2((X_addback[k]),(Y_addback[k]),(Z_addback[k]));
	            TVector3 V3(0,0,390);
	            Float_t addback_gAngle_Origin = (V3.Angle(V2))*180./M_PI;
	
	            //This assumes that recoil hits the (icPositionX,icPositionY,390) in IC OR (Recoil_X,Recoil_Y,390) from calculation (CHOOSE ONE)
	            //TVector3 V4(icPositionX-TargetX, icPositionY-TargetY, 390-TargetZ); //This uses (icPositionX,icPositionY,390)
                TVector3 V4(Recoil_X-TargetX, Recoil_Y-TargetY, 390-TargetZ); //This uses (Recoil_X,Recoil_Y,390)
	            Float_t addback_gAngle_IC = (V4.Angle(V1))*180./M_PI;
                
	            //Doppler Correcting here
                Float_t addback_dGamma_Origin = dop(addback_gAngle_Origin,beta,xtals_addback[k]);
                Float_t addback_dGamma_IC = dop(addback_gAngle_IC,beta,xtals_addback[k]);
                

                //Using the (d,p) events cut in kinematics
                //if(!kinCut->IsInside(Q5Angle, (QQQ5RingEnergyTarget*1000))) continue;

                //Fill your histogram after addback
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && timeStamp>0 && GRETINATimeStamp>0 && foundGRETINA && (timeStamp-GRETINATimeStamp)>20. && (timeStamp-GRETINATimeStamp)<60.) dopAddGenVQen_Timecut->Fill(addback_dGamma_Origin,(RingExcitationEnergy*1000.));     
                if (QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && timeStamp>0 && GRETINATimeStamp>0 && foundGRETINA && (timeStamp-GRETINATimeStamp)>20. && (timeStamp-GRETINATimeStamp)<60.) dopAddGenVQen_Timecut_Q5->Fill(addback_dGamma_Origin,(RingExcitationEnergy*1000.));     

                //Fill your particle histogram with cut on gamma rays
                if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && addback_dGamma_Origin>=1250. && addback_dGamma_Origin<=1290.) Q5Ex_TDC_1266->Fill(RingExcitationEnergy);
                if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && addback_dGamma_Origin>=3130. && addback_dGamma_Origin<=3150.) Q5Ex_TDC_3134->Fill(RingExcitationEnergy);
                if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && addback_dGamma_Origin>=2010. && addback_dGamma_Origin<=2040.) Q5Ex_TDC_2028->Fill(RingExcitationEnergy);
                if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && addback_dGamma_Origin>=2120. && addback_dGamma_Origin<=2190.) Q5Ex_TDC_2148->Fill(RingExcitationEnergy);
                if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100. && addback_dGamma_Origin>=3486. && addback_dGamma_Origin<=3526.) Q5Ex_TDC_3506->Fill(RingExcitationEnergy);
            } 

//---------------------------------------------------------End of Addback----------------------------------------------------------------------------//
	        
            //Fill your QQQ5 Histograms Here
            
	        if(QQQ5Det[j]==1 && QQQ5Sector[j]==3)continue;

	        // QQQ5 Stuffs
	        if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) Q5Ex_TDC->Fill(RingExcitationEnergy);
            if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) Q5EvA_TDC->Fill(Q5Angle,(1000.*QQQ5RingEnergyTarget));
	        if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.) Q5_RvS_TDC->Fill(QQQ5RingEnergy[j]/1000.,QQQ5SectorEnergy[j]/1000.);

            //Using More Cuts in QQQ5
            //Front vs Back cut in QQQ5 Detector
            if(!QQQ5Cut_FvB->IsInside(QQQ5RingEnergy[j]/1000.,QQQ5SectorEnergy[j]/1000.)) continue;
            if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.)Q5_RvS_TDC_cut->Fill(QQQ5RingEnergy[j]/1000.,QQQ5SectorEnergy[j]/1000.);

            //(d,p) events cut in kinematics
            if(!kinCut->IsInside(Q5Angle, (QQQ5RingEnergyTarget*1000))) continue;
            if(QQQ5Upstream[j] && tdcIC>700. && tdcIC<1100.)Q5EvA_TDC_kinCut->Fill(Q5Angle,(1000.*QQQ5RingEnergyTarget));
        }//End of QQQ5 Multiplicity

//--------------------------------------------------------------QQQ5 Multiplicity Loop ends---------------------------------------------------------------------------------------------//

//---------------------------------------------------------------SX3 Multiplicity Loop Start--------------------------------------------------------------------------------------------//

        //Loop over the SX3 Multiplicity
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
            //Float_t ScaleX = 0.3; 
			//Float_t ScaleY = 0.36; 
            Float_t ScaleX = 0; 
			Float_t ScaleY = 0;
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
            SectorExcitationEnergy = SectorExcitationEnergy * 0.980736 + (-0.32265);

            //Finding the beta for doppler correction
            float QVal = Qgs - SectorExcitationEnergy;
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

//-----------------------------------------------------------------Crystal Multiplicity Loop Start--------------------------------------------------------------------------------------------------------//

            //Xtal Singles
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
                    X_xtals[k] = X_xtals[k] - GretinaOffset_X1; // Beam Left Shift
                } else{
                    X_xtals[k] = X_xtals[k] + GretinaOffset_X2; //Beam Right Shift
                }
                Y_xtals[k] = Y_xtals[k] - GretinaOffset_Y; //ORRUBA  higher than GRETINA
                Z_xtals[k] = Z_xtals[k] - GretinaOffset_Z; //ORRUBA downstream of GRETINA

	            //Finding the gammaAngle for doppler correction
	            //If we assume pin-point beam at origin
                TVector3 V1((X_xtals[k]-TargetX),(Y_xtals[k]-TargetY),(Z_xtals[k]-TargetZ)); 
	            TVector3 V2((X_xtals[k]),(Y_xtals[k]),(Z_xtals[k]));
	            TVector3 V3(0,0,390);
	            Float_t gAngle_Origin = (V3.Angle(V2))*180./M_PI;
                
                 //This assumes that recoil hits the (icPositionX,icPositionY,390) in IC OR (Recoil_X,Recoil_Y,390) from calculation (CHOOSE ONE)
	            TVector3 V4(icPositionX-TargetX, icPositionY-TargetY, 390-TargetZ); //This uses (icPositionX,icPositionY,390)
                //TVector3 V4(Recoil_X-TargetX, Recoil_Y-TargetY, 390-TargetZ); //This uses (Recoil_X,Recoil_Y,390)
	            Float_t gAngle_IC = (V4.Angle(V1))*180./M_PI;
                
	            //Doppler Correcting here
                Float_t dGamma_Origin = dop(gAngle_Origin,beta,xtals_cc[k]);
                Float_t dGamma_IC = dop(gAngle_IC,beta,xtals_cc[k]);

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
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950.) GenVQen_SX3->Fill(xtals_cc[k],(SectorExcitationEnergy*1000.));

                //Using Other cuts here
                //Front vs Back cut in SX3 Detectors
                if(!SX3Cut_FvB->IsInside(Energy,SectorEnergy)) continue;
                //(d,p) kinematics cut 
                //if(!kinCut->IsInside(SX3Angle, SectorEnergyTarget)) continue;
                
                //Filling the histogram after the Front Vs Back cut and the (d,p) kinematics cut
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950. && timeStamp>0 && GRETINATimeStamp>0 && foundGRETINA && (timeStamp-GRETINATimeStamp)>20. && (timeStamp-GRETINATimeStamp)<60.) dopGenVQen_Timecut->Fill(dGamma_Origin,(SectorExcitationEnergy*1000.));
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950. && timeStamp>0 && GRETINATimeStamp>0 && foundGRETINA && (timeStamp-GRETINATimeStamp)>20. && (timeStamp-GRETINATimeStamp)<60.) dopGenVQen_Timecut_SX3->Fill(dGamma_Origin,(SectorExcitationEnergy*1000.));
            }//End of XTAL Singles    

			
        //---------------------------------------------------------- XTAL ADDBACK -------------------------------------------------------------------------------------//	       
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
            }//End of loop over k   

//--------------------------------------------------------------XTAL ADDBACK ENDS & DOPPLER CORRECTING "ADDBACK-ED" GAMMAS STARTS HERE--------------------------------------------------------------------------------------------------//
            Float_t addback_gAngle_Origin[gmult];
            for(Int_t k=0; k<gmult; k++){
                //If we assume pin-point beam at origin
                TVector3 V1((X_addback[k]-TargetX),(Y_addback[k]-TargetY),(Z_addback[k]-TargetZ)); 
	            TVector3 V2((X_addback[k]),(Y_addback[k]),(Z_addback[k]));
	            TVector3 V3(0,0,390);
	            Float_t addback_gAngle_Origin = (V3.Angle(V2))*180./M_PI;
	
	            //This assumes that recoil hits the (icPositionX,icPositionY,390) in IC OR (Recoil_X,Recoil_Y,390) from calculation (CHOOSE ONE)
	            //TVector3 V4(icPositionX-TargetX, icPositionY-TargetY, 390-TargetZ); //This uses (icPositionX,icPositionY,390)
                TVector3 V4(Recoil_X-TargetX, Recoil_Y-TargetY, 390-TargetZ); //This uses (Recoil_X,Recoil_Y,390)
	            Float_t addback_gAngle_IC = (V4.Angle(V1))*180./M_PI;
                
	            //Doppler Correcting here
                Float_t addback_dGamma_Origin = dop(addback_gAngle_Origin,beta,xtals_addback[k]);
                Float_t addback_dGamma_IC = dop(addback_gAngle_IC,beta,xtals_addback[k]);
                

                //Using the (d,p) events cut in kinematics
                //if(!kinCut->IsInside(Q5Angle, (QQQ5RingEnergyTarget*1000))) continue;

                //Fill your histogram after addback
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950. && timeStamp>0 && GRETINATimeStamp>0 && foundGRETINA && (timeStamp-GRETINATimeStamp)>20. && (timeStamp-GRETINATimeStamp)<60.) dopAddGenVQen_Timecut->Fill(addback_dGamma_Origin,(SectorExcitationEnergy*1000.));
                if (SX3Upstream[j] && tdcIC>600. && tdcIC<950. && timeStamp>0 && GRETINATimeStamp>0 && foundGRETINA && (timeStamp-GRETINATimeStamp)>20. && (timeStamp-GRETINATimeStamp)<60.) dopAddGenVQen_Timecut_SX3->Fill(addback_dGamma_Origin,(SectorExcitationEnergy*1000.));
                
                //Fill your particle histogram with cut on gamma rays
			    if(SX3Upstream[j] && tdcIC>600. && tdcIC<950. && addback_dGamma_Origin>=1250. && addback_dGamma_Origin<=1290.)SX3Ex_TDC_1266->Fill(SectorExcitationEnergy);
 			    if(SX3Upstream[j] && tdcIC>600. && tdcIC<950. && addback_dGamma_Origin>=3110. && addback_dGamma_Origin<=3140.)SX3Ex_TDC_3134->Fill(SectorExcitationEnergy);
			    if(SX3Upstream[j] && tdcIC>600. && tdcIC<950. && addback_dGamma_Origin>=2010. && addback_dGamma_Origin<=2045.)SX3Ex_TDC_2028->Fill(SectorExcitationEnergy);
			    if(SX3Upstream[j] && tdcIC>600. && tdcIC<950. && addback_dGamma_Origin>=2135. && addback_dGamma_Origin<=2165.)SX3Ex_TDC_2148->Fill(SectorExcitationEnergy);
			    if(SX3Upstream[j] && tdcIC>600. && tdcIC<950. && addback_dGamma_Origin>=3486. && addback_dGamma_Origin<=3526.)SX3Ex_TDC_3506->Fill(SectorExcitationEnergy);
            } 

//---------------------------------------------------------End of Addback----------------------------------------------------------------------------//
            //Fill your SX3 histograms here
			if(SX3Upstream[j] && tdcIC>600. && tdcIC<950.)Q5EvA_TDC->Fill(SX3Angle,SectorEnergyTarget);
			if(SX3Upstream[j] && tdcIC>600. && tdcIC<950.)SX3Ex_TDC->Fill(SectorExcitationEnergy);
            if(SX3Upstream[j] && tdcIC>600. && tdcIC<950.)SX3_FvB_TDC->Fill(Energy,SectorEnergy);
			
            //Using More Cuts in QQQ5
            //Front vs Back cut in SX3 Detectors
            if(!SX3Cut_FvB->IsInside(Energy,SectorEnergy)) continue;
            if(SX3Upstream[j] && tdcIC>600. && tdcIC<950.)SX3_FvB_TDC_cut->Fill(Energy,SectorEnergy);

            //Using the (d,p) kinematics cut 
            if(!kinCut->IsInside(SX3Angle, SectorEnergyTarget)) continue;
            if(SX3Upstream[j] && tdcIC>600. && tdcIC<950.)Q5EvA_TDC_kinCut->Fill(SX3Angle,SectorEnergyTarget);
        }//End of Loop over SX3 Multiplicity    
    }// End of event by event analysis

//----------------------------------------------------Event by event analysis ends here--------------------------------------------------------------------------------------------//
  
    outputFile->cd();

    // Write your histograms here
    //IC stuffs
    hICdEE->Write();
    hICdEE_cut->Write();
    hICdEE_kincut->Write();

    //QQQ5 Stuffs
    Q5EvA_TDC->Write();
    Q5Ex_TDC->Write();
    Q5EvA_TDC_kinCut->Write();
    Q5_RvS_TDC->Write();
    Q5_RvS_TDC_cut->Write();

    //SX3 Stuffs
    SX3_FvB_TDC->Write();
    SX3_FvB_TDC_cut->Write();
    SX3Ex_TDC->Write();

    //Gamma Stuffs
    GenVQen->Write();
    GenVQen_Q5->Write();
    GenVQen_SX3->Write();
    dopGenVQen_Timecut->Write();
    dopGenVQen_Timecut_Q5->Write();
    dopGenVQen_Timecut_SX3->Write();
   
    //Addback Stuffs
    dopAddGenVQen_Timecut->Write();
    dopAddGenVQen_Timecut_Q5->Write();
    dopAddGenVQen_Timecut_SX3->Write();

    //Particle Stuffs With Gamma cuts
    Q5Ex_TDC_1266->Write();    
    SX3Ex_TDC_1266->Write();
    Q5Ex_TDC_3134->Write();    
    SX3Ex_TDC_3134->Write();
    Q5Ex_TDC_2028->Write();    
    SX3Ex_TDC_2028->Write();
    Q5Ex_TDC_2148->Write();    
    SX3Ex_TDC_2148->Write();
    Q5Ex_TDC_3506->Write();    
    SX3Ex_TDC_3506->Write();

    outputFile->Close();

}

//-------------------------------------------------------------End of the analysis loop-----------------------------------------------------------------------------------------------//
