#include "run.hh"
#include "ions.hh"

MyRunAction::MyRunAction(): fRun(nullptr), fNumOfEvent(0), fRunID(0)
{

	minlogE = 1.0;
	maxlogE = 5.0;
 	NbinsE = 100;
	
	
	// Get the current time in seconds since the epoch
	std::time_t currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

	// Convert the time to a struct tm in the local timezone
	std::tm* timeInfo = std::localtime(&currentTime);

	// Create a stringstream to format the date and time
	std::stringstream ss;

	// Format the date and time as "YYYYMMDD-hhmmss"
	ss << std::put_time(timeInfo, "%Y%m%d-%H%M%S");

        std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
        auto duration = now.time_since_epoch();
        auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);

	//G4String rd_label = G4String(std::to_string(CLHEP::RandFlat::shootInt(1e9)));
	G4String rd_label = G4String(std::to_string(nanoseconds.count()%100000000));

	labelCSV = G4String (ss.str()) + G4String("-") + rd_label;

	//labelCSV = G4String(ss.str());
	DefineCommands();
}

MyRunAction::~MyRunAction()
{
	delete fMessenger;
}


G4Run* MyRunAction::GenerateRun()
{
	const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	
	phantomType = detectorConstruction->GetPhantomType();
	//SetResultsDirectory("/results/"+phantomType+G4String("/"));
	G4cout << "GenerateRun " << phantomType << G4endl;
	if ((phantomType=="ICRP145") || (phantomType=="BDRTOG4")) 
	{	
		fRun = new TETRun();
		return fRun;
	}
	else
	{
		iRun = new ISRun();
		return iRun;
	}
}


void MyRunAction::UpdateFlux(G4double ikE, G4double secondarykE, G4String secondaryName)
{
  	G4int energyBin = G4int((log10(ikE)-minlogE)/(maxlogE-minlogE)*NbinsE);
	flux[energyBin].Update(secondarykE, secondaryName);
}

void MyRunAction::IterIncident(G4double ikE)
{
  	G4int energyBin = G4int((log10(ikE)-minlogE)/(maxlogE-minlogE)*NbinsE);
	
	if (flux.count(energyBin)==0)
	{
		//G4String csvSource = this->GetCSVFluxName(primaryName);
		flux[energyBin].FirstUse();
	}
	flux[energyBin].IterIncident();
}


void MyRunAction::BeginOfRunAction(const G4Run* aRun)
{
	G4cout << "Begin run action" << G4endl;
	auto man = G4AnalysisManager::Instance();
	
	const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	
	//phantomType = detectorConstruction->GetPhantomType();

	if ((phantomType=="ICRP145") || (phantomType=="BDRTOG4"))
	{
		fTetData = detectorConstruction->GetICRP145phantom();
		auto massMap = fTetData->GetMassMap();
	}
	const MyPrimaryGenerator *generator = static_cast<const MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	//std::vector<G4String> particleList;

	//for (G4int i=0;i<generator->GetNIons();i++) 
	//{
		//G4cout << i << " " << generator->GetIonName(i) << G4endl;
	//	particleList.push_back(generator->GetIonName(i));
	//}
	G4String primaryParticle = generator->GetIonName();


	man->CreateNtuple("Doses","Doses");
	man->CreateNtupleIColumn("organId");
	man->CreateNtupleIColumn("eBin");
	//man->CreateNtupleIColumn(G4String("N"));
	man->CreateNtupleDColumn(G4String("DE"));
	man->CreateNtupleDColumn(G4String("AD"));
	man->FinishNtuple(0);


        man->CreateNtuple("InnerFlux","InnerFlux");
        man->CreateNtupleIColumn("ikE");
        man->CreateNtupleSColumn("Oparticle");
        man->CreateNtupleDColumn("okE");
        man->CreateNtupleDColumn("count");
        man->FinishNtuple(1);

	man->CreateNtuple("N","N");
	man->CreateNtupleIColumn("ikE");
	man->CreateNtupleIColumn("N");
	man->FinishNtuple(2);
	
	if ((phantomType=="ICRP145") || (phantomType=="BDRTOG4")) {
		// print the progress at the interval of 10%
		fNumOfEvent=aRun->GetNumberOfEventToBeProcessed();
		G4RunManager::GetRunManager()->SetPrintProgress(int(fNumOfEvent*0.1));
	}
	
	G4cout << "Run created" << G4endl;
}

void MyRunAction::EndOfRunAction(const G4Run* aRun)
{


	const MyPrimaryGenerator *generator = static_cast<const MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	
	auto man = G4AnalysisManager::Instance();

	for (int i=0;i<100;i++)
	{
		man->FillNtupleIColumn(2,0,i);
		man->FillNtupleIColumn(2,1,generator->GetNGenerated(i));
		man->AddNtupleRow(2);
	}

	std::map<G4int, TotalFlux>::iterator itF;

	std::map<G4String,Ion> Ions = getIons();

        // Record fluxes
        for (itF = flux.begin(); itF != flux.end(); itF++)
        {
                G4int ikE = itF->first;
                G4String oParticle;
                //G4int iZ = Ions[iParticle].getZ();
		G4int iZ = 1;

                //man->FillNtupleSColumn(4,0,iParticle);
                //man->FillNtupleIColumn(4,1,iZ);
                //man->FillNtupleDColumn(4,2,generator->GetTotalParticleNumber(iZ-1));
                //man->AddNtupleRow(4);

                std::map<G4String,ParticleSpectra> ofluxes = itF->second.GetFluxes();
                std::map<G4String,ParticleSpectra>::iterator itoF;

                for (itoF = ofluxes.begin(); itoF != ofluxes.end();itoF++)
                {
                        oParticle = itoF->first;
                        for (G4int oebin=0; oebin<itoF->second.GetNbin();oebin++)
                        {
                                man->FillNtupleIColumn(1,0,ikE);
                                man->FillNtupleSColumn(1,1,oParticle);
                                man->FillNtupleDColumn(1,2,oebin);
                                man->FillNtupleDColumn(1,3,itoF->second.GetBin(oebin));
                                man->AddNtupleRow(1);
                        }
                }

        }

	man->Write();
	man->CloseFile();

 	// print the result only in the Master
 	if(!isMaster) return;

 	// get the run ID
 	fRunID = aRun->GetRunID();
	
}

void MyRunAction::SetResultsDirectory(G4String dir)
{
	const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	
	phantomType = detectorConstruction->GetPhantomType();
	auto man = G4AnalysisManager::Instance();
	G4cout << "Set Dir Phantom " << phantomType << " " << dir << G4endl;
	//result_dir = "../results/"+phantomType+G4String("/")+dir+G4String("/");
	result_dir = dir+G4String("/");
	//G4cout << result_dir << G4endl;
	
	createDirectories(result_dir);
	
	//auto man = G4AnalysisManager::Instance();
	man->SetDefaultFileType("csv");
	//man->SetVerboseLevel(1);
	//G4cout << result_dir+labelCSV+ G4String(".csv") << G4endl;
	man->OpenFile(result_dir+labelCSV+ G4String(".csv"));
	//G4cout << result_dir+labelCSV+ G4String(".csv") << G4endl;
	//man->OpenFile(G4String("output.csv"));
}

void MyRunAction::DefineCommands()
{
	fMessenger = new G4GenericMessenger(this, "/SIM/scoring/","Set scoring folder");

	auto& resdirCmd = fMessenger->DeclareMethod("resDir",&MyRunAction::SetResultsDirectory,"Folder for results");
	resdirCmd.SetParameterName("resDir", true);
	resdirCmd.SetDefaultValue("../results/"+phantomType+G4String("/"));
}


void MyRunAction::PrintResult(std::ostream &out)
{

if (phantomType!="ICRP145") {return;}
 // Print run result
 //
 using namespace std;
 EMAP edepMap = *fRun->GetEDEMap();

    out << G4endl
	 << "=====================================================================" << G4endl
	 << " Run #" << fRunID << " / Number of event processed : "<< fNumOfEvent    << G4endl
	 << "=====================================================================" << G4endl
	 << "organ ID| "
	 << setw(19) << "Organ Mass (g)"
         << setw(19) << "Dose (Gy/source)"
	 //<< setw(19) << "Relative Error" 
	 << setw(34) << "Name of organ" << G4endl;

    out.precision(3);
    auto massMap = fTetData->GetMassMap();
    for(auto itr : massMap){
    		//G4String id_map = G4String(std::to_string(itr.first));
    		
		G4double meanDose    = edepMap[0][itr.first]  / itr.second / fNumOfEvent;
		//G4double squareDose =  (edepMap[itr.first].second)/ (itr.second*itr.second);
		//G4double variance    = ((squareDose/fNumOfEvent) - (meanDose*meanDose))/fNumOfEvent;
		G4double relativeE(1);
		//if(meanDose > 0) relativeE   = sqrt(variance)/meanDose;

		out << setw(8)  << itr.first << "| "
			<< setw(19) << fixed      << itr.second/g;
		out	<< setw(19) << scientific << meanDose/(joule/kg);
		//out	<< setw(19) << fixed      << relativeE ; 
		out	<< setw(34) << fTetData->GetMaterial(itr.first)->GetName() << G4endl;
	}
	out << "=====================================================================" << G4endl << G4endl;
}

void createDirectories(const std::string& path) {
    std::string currentPath = "";
    for (char c : path) {
        if (c != '/') {
            currentPath += c;
        } else {
            currentPath += '/';
            if (!std::filesystem::exists(currentPath)) {
                std::filesystem::create_directory(currentPath);
            }
        }
    }
    if (!std::filesystem::exists(currentPath)) {
        std::filesystem::create_directory(currentPath);
    }
}
