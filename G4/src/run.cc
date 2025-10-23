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


        /*man->CreateNtuple("InnerFlux","InnerFlux");
        man->CreateNtupleIColumn("ikE");
        man->CreateNtupleSColumn("Oparticle");
        man->CreateNtupleDColumn("okE");
        man->CreateNtupleDColumn("count");
        man->FinishNtuple(1);*/

	man->CreateNtuple("N","N");
	man->CreateNtupleIColumn("ikE");
	man->CreateNtupleIColumn("N");
	man->FinishNtuple(1);
	
	if ((phantomType=="ICRP145") || (phantomType=="BDRTOG4")) {
		// print the progress at the interval of 10%
		fNumOfEvent=aRun->GetNumberOfEventToBeProcessed();
		G4RunManager::GetRunManager()->SetPrintProgress(int(fNumOfEvent*0.1));
	}
	
}

void MyRunAction::EndOfRunAction(const G4Run* aRun)
{


	G4cout << "EndOfRunAction" << G4endl;
	const MyPrimaryGenerator *generator = static_cast<const MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	
	auto man = G4AnalysisManager::Instance();


	//std::map<G4int,G4int> Npart = generator->GetNGenerated();
	for (const auto &[key, value]: generator->GetNGenerated())
	{
		man->FillNtupleIColumn(1,0,key);
		man->FillNtupleIColumn(1,1,value);
		man->AddNtupleRow(1);
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
