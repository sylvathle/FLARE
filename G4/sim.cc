#include <iostream>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "FTFP_BERT_HP.hh"
#include "QGSP_INCLXX_HP.hh"
#include "QBBC.hh"
#include "Shielding.hh"
#include "QGSP_BERT.hh"

#include "geometry.hh"
//#include "sphere_geometry.hh"
//#include "physics.hh"
#include "action.hh"

#include "shieldmaterial.hh"

int main(int argc, char** argv)
{
	time_t start,end;
	time(&start);

	G4RunManager *runManager = new G4RunManager();

	runManager->SetVerboseLevel(0);
	G4UIExecutive *ui = nullptr;

	if (argc<=1) {ui = new G4UIExecutive(argc,argv);}


	// Decide Physics list to use.
	runManager->SetUserInitialization(new Shielding(0));
	//runManager->SetUserInitialization(new FTFP_BERT_HP(0));
	//runManager->SetUserInitialization(new QGSP_INCLXX_HP(0));
	//runManager->SetUserInitialization(new QBBC(0));
	

	runManager->SetUserInitialization(new MyGeometry(ui));
	runManager->SetUserInitialization(new MyActionInitialization());
	
	

	//const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (runManager->GetUserDetectorConstruction());
	//G4String phantomType = detectorConstruction->GetPhantomType();
	
	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	

	if ( ! ui)
	{
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);

	}
	else
	{
		G4VisManager *visManager = new G4VisExecutive;
		visManager->SetVerboseLevel(0);
		visManager->Initialize();
		UImanager->ApplyCommand(G4String("/control/execute vis.mac"));
		ui->SessionStart();
		delete visManager;
	}

	delete runManager;
	delete ui;

	time(&end);

	double time_taken = double(end-start);
	G4cout << "Execution time: " << time_taken << std::setprecision(5) << "sec" << G4endl;
}
