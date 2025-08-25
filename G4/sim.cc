#include <iostream>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//#include "QBBC_HP.hh"
#include "QBBC.hh"
#include "Shielding.hh"
#include "QGSP_BERT.hh"

#include "geometry.hh"
//#include "sphere_geometry.hh"
#include "physics.hh"
#include "action.hh"

#include "shieldmaterial.hh"

int main(int argc, char** argv)
{
	time_t start,end;

	time(&start);

	//#ifndef G4MULTITHREADED
		G4RunManager *runManager = new G4RunManager();
	
	//#else
	//	G4MTRunManager *runManager = new G4MTRunManager();
	//#endif

	runManager->SetVerboseLevel(0);
	G4UIExecutive *ui = nullptr;

	//G4String *shieldmaterial= new G4String("G4_Galactic");
	//G4String *shieldmaterialname= new G4String("G4_Al");
	//ShieldMaterial* shieldmaterial = new ShieldMaterial(shieldmaterialname);
	//G4cout << "argc " << argc << G4endl;
	//if (argc==3)
	//{
	//	
	//	shieldmaterialname= new G4String(argv[2]);
	//	shieldmaterial=new ShieldMaterial(shieldmaterialname);
	//	G4cout << shieldmaterial << G4endl;
	//}

	if (argc<=1) {ui = new G4UIExecutive(argc,argv);}

	//G4String phantomType = "ICRP145"; //"ICRP145" "IcruSphere"
	//phantomType = "IcruSphere"; //"ICRP145" "IcruSphere"

	runManager->SetUserInitialization(new Shielding(0));
	runManager->SetUserInitialization(new MyGeometry(ui));
	runManager->SetUserInitialization(new MyActionInitialization());
	
	//runManager->SetUserInitialization(new MyPhysicsList());
	//runManager->SetUserInitialization(new QBBC(0));
	
	

	//G4VModularPhysicsList* physics = new QGSP_BERT();
	//G4VModularPhysicsList* physics = new QBBC();
    	//physics->RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP());
	//runManager->SetUserInitialization(physics);
	

	const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (runManager->GetUserDetectorConstruction());
	G4String phantomType = detectorConstruction->GetPhantomType();
	
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
		//G4cout << G4String("/control/execute vis")+phantomType+G4String(".mac") << G4endl;
		//UImanager->ApplyCommand(G4String("/control/execute vis")+phantomType+G4String(".mac"));
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
