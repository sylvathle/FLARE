#include "geometry.hh"
//#include "CADMesh.hh"

//MyGeometry::MyGeometry(): phantomType("IcruSphere"), innerRadius(0.2*m), thickVest(0.1*m)
MyGeometry::MyGeometry(): phantomType("IcruSphere"), thickVest(0.1*m)
//MyGeometry::MyGeometry(): regoWidth(0.05*m), aluWidth(0.01*m), rPhantom(0.15*m),phantomType("ICRP145")
{
	materials = new Materials();
	DefineCommands();
}

MyGeometry::MyGeometry(G4UIExecutive *ui_): rPhantom(0.15*m),phantomType("IcruSphere"), thickVest(0.1*m)
//MyGeometry::MyGeometry(G4UIExecutive *ui_): rPhantom(0.15*m),phantomType("IcruSphere"), innerRadius(0.2*m), thickVest(0.1*m)
//MyGeometry::MyGeometry(G4UIExecutive *ui_): regoWidth(0.05*m), aluWidth(0.01*m), rPhantom(0.15*m),phantomType("ICRP145")
{

	ui = ui;
	materials = new Materials();
	DefineCommands();
}
	

MyGeometry::~MyGeometry()
{
}	

G4VPhysicalVolume *MyGeometry::Construct()
{

	// Import TET model for the ICRP 145 HP
	if (phantomType == "ICRP145")
	{

		fTetData = new TETModelImport(true,true,ui);
		fPhantomSize     = fTetData -> GetPhantomSize();
		fPhantomBoxMin   = fTetData -> GetPhantomBoxMin();
		fPhantomBoxMax   = fTetData -> GetPhantomBoxMax();
		fNOfTetrahedrons = fTetData -> GetNumTetrahedron();
	}

	if (phantomType == "BDRTOG4")
	{
		fTetData = new TETModelImport(true,ui,csvBodies);
		fPhantomSize     = fTetData -> GetPhantomSize();
		fPhantomBoxMin   = fTetData -> GetPhantomBoxMin();
		fPhantomBoxMax   = fTetData -> GetPhantomBoxMax();
		fNOfTetrahedrons = fTetData -> GetNumTetrahedron();
	}
	
	// Dimensions of the world
	G4double xWorld = 10.0*m;
	G4double yWorld = 10.0*m;
	G4double zWorld = 10.0*m;

	// Create world of simulation.
	solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
	logicWorld = new G4LogicalVolume(solidWorld,materials->GetMaterial("G4_Galactic"),"logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, true);
	
	heightModule = 1.0*m;
	//radiusModule = 3.346/6.0*m;
	
	radiusModule = 1.2*m;
	thickModule = 4*mm;
	
	//solidAir =  new G4Tubs("solidAir",0.0*m,radiusModule,heightModule,0*deg,360*deg);
	solidAir = new G4Sphere("solidDome3", 0.,  radiusModule, 0*deg, 360*deg, 0*deg, 180*deg);
	logicAir = new G4LogicalVolume(solidAir,materials->GetMaterial("G4_Galactic"),"logicAir");
	physAir = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicAir, "physAir", logicWorld, false, 0, true);
	
        //innerRadius = fPhantomSize.y()-1.0*cm;
        //innerRadius = 1.0*m
        //out_radius1 = innerRadius+thickVest;
        G4double heightVest = 0.36*m;
        //G4ThreeVector placementCylinder(0., 0.05*m, 0.28*m);
        G4double cylinderHeight = 0.14*m;
        G4ThreeVector placementCylinder(0., 0.*m, 0.);
	G4RotationMatrix* rot = new G4RotationMatrix();

	if (phantomType=="ICRP145")
	{

		shift = G4ThreeVector(0.0*m,0.0*m,0.0*m);
		
		//G4cout << "Geometry, setting up ICRP145" << G4endl; 
 		containerSolid = new G4Ellipsoid("phantomBox", fPhantomSize.x()/2 + 4.5*cm,
					           fPhantomSize.y()/2 + 7.5*cm,
						   fPhantomSize.z()/2 + 23.*cm, -fPhantomSize.z()/2+1.0*cm,fPhantomSize.z()/2+2.0*cm);


 		fContainer_logic = new G4LogicalVolume(containerSolid, materials->GetMaterial("G4_AIR"), "logicPhantom");
	
		new G4PVPlacement(rot, G4ThreeVector(0.0*m,0.0*m,0.0*m), fContainer_logic, "PhantomPhysical", logicAir, false, 0, true);		

		fContainer_logic->SetOptimisation(TRUE);
		fContainer_logic->SetSmartless( 0.5 ); // for optimization (default=2)

		// Define the tetrahedral mesh phantom as a parameterised geometry
		//
		// solid and logical volume to be used for parameterised geometry
		G4VSolid* tetraSolid = new G4Tet("TetSolid", shift,  G4ThreeVector(1.*cm,0,0), G4ThreeVector(0,1.*cm,0),  G4ThreeVector(0,0,1.*cm));
	
		fTetLogic = new G4LogicalVolume(tetraSolid, materials->GetMaterial("G4_Galactic"), "TetLogic");
		
		// physical volume (phantom) constructed as parameterised geometry
		new G4PVParameterised("wholePhantom",fTetLogic,fContainer_logic, kUndefined, fTetData->GetNumTetrahedron(), new TETParameterisation(fTetData));

		
	}
	else if (phantomType=="BDRTOG4")
	{
		shift = G4ThreeVector(0.0*m,0.0*m,0.0*m);

		// Define the tetrahedral mesh phantom as a parameterised geometry
		//
		// solid and logical volume to be used for parameterised geometry
		G4VSolid* tetraSolid = new G4Tet("TetSolid", shift,  G4ThreeVector(1.*cm,0,0), G4ThreeVector(0,1.*cm,0),  G4ThreeVector(0,0,1.*cm));	
		fTetLogic = new G4LogicalVolume(tetraSolid, materials->GetMaterial("G4_Galactic"), "TetLogic");
		
		// physical volume (phantom) constructed as parameterised geometry
		new G4PVParameterised("wholePhantom",fTetLogic,logicAir, kUndefined, fTetData->GetNumTetrahedron(), new TETParameterisation(fTetData));
		
	}
	else if (phantomType=="IcruSphere")
	{
		solidPhantom = new G4Sphere("solidPhantom", 0., rPhantom, 0, 360*deg, 0, 180*deg);
		logicPhantom = new G4LogicalVolume(solidPhantom, materials->GetMaterial("IcruMat"), "logicPhantom");
		physPhantom = new G4PVPlacement(0, G4ThreeVector(0.,0.*m,0.), logicPhantom, "physPhantom", logicWorld, false, 0, true);
	}
	
	return physWorld;
}

void MyGeometry::ConstructSDandField()
{
	G4SDManager* pSDman = G4SDManager::GetSDMpointer();
	G4String phantomSDname = "PhantomSD";
	// MultiFunctional detector
 	//auto* MFDet = new MyMultiDetector(phantomSDname);
 	auto* MFDet = new G4MultiFunctionalDetector(phantomSDname);
 	pSDman->AddNewDetector( MFDet );
		

	if ((phantomType=="ICRP145") || (phantomType=="BDRTOG4"))
	{
		// scorer for energy depositon in each organ
		MFDet->RegisterPrimitive(new TETPSEnergyDeposit("EDE", fTetData));
		//MFDet->RegisterPrimitive(new TETPSEnergyDeposit("Dose", fTetData));

	//	// attach the detector to logical volume for parameterised geometry (phantom geometry)
		SetSensitiveDetector(fTetLogic, MFDet);
		
	}
	else
	{
		//G4cout << "Sensitive dec " << phantomType << G4endl;
		MFDet->RegisterPrimitive(new ISPSEnergyDeposit("EDE"));

		SetSensitiveDetector(logicPhantom, MFDet);
	}

	G4cout << "ConstructSDandField" << G4endl;

	G4MultiFunctionalDetector* myScorer = new G4MultiFunctionalDetector("myCellScorer");
	G4SDManager::GetSDMpointer()->AddNewDetector(myScorer);
	G4VPrimitiveScorer* totalSurfFlux = new SBG4PSSphereSurfaceFlux("TotalSurfFlux",0);
	myScorer->RegisterPrimitive(totalSurfFlux);
	logicAir->SetSensitiveDetector(myScorer);

}

void MyGeometry::SetHumanPhantom(G4String phantom_) {phantomType = phantom_;}

void MyGeometry::SetModule(G4double thickness=4*mm)
{
	thickModule = thickness;
	if (thickness>0.1*mm)
	{
		G4cout << "Module thickness = " << thickness << G4endl;
		solidModule = new G4Sphere("solidDome3", radiusModule, radiusModule+thickness, 0*deg, 360*deg, 0*deg, 180*deg);
		logicModule = new G4LogicalVolume(solidModule,materials->GetMaterial("G4_Al"),"logicModule");
		physModule = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicModule, "physModule", logicWorld, false, 0, true);
		
		//static_cast<MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())->SetBeamRadius(thickness+rsource+0.2*m);
		//generator->SetBeamRadius(thickness+rsource+0.2*m);
			
		/*solidModuleUp = new G4Tubs("solidModuleUp", 0.0*m, radiusModule+thickModule, thickModule/2.0, 0*deg, 360*deg);
		logicModuleUp = new G4LogicalVolume(solidModuleUp,materials->GetMaterial("G4_Al"),"logicModuleUp");
		physModuleUp = new G4PVPlacement(0, G4ThreeVector(0.,0.,heightModule+thickModule/2.0), logicModuleUp, "physModuleUp", logicWorld, false, 0, true);

		solidModuleDown = new G4Tubs("solidModuleDown", 0.0*m, radiusModule+thickModule, thickModule/2.0, 0*deg, 360*deg);
		logicModuleDown = new G4LogicalVolume(solidModuleDown,materials->GetMaterial("G4_Al"),"logicModuleDown");
		physModuleDown = new G4PVPlacement(0, G4ThreeVector(0.,0.,-heightModule-thickModule/2.0), logicModuleUp, "physModuleDown", logicWorld, false, 0, true);*/
	
		logicAir->SetMaterial(materials->GetMaterial("G4_AIR"));
	}
}

void MyGeometry::SetThickVest(G4double t1)
{	
	if (t1==0)
	{
		G4double thick = 1*mm;
		delete solidVest,logicVest,physVest,solidFilledVest;				
	}

	else
	{
		G4double thick = t1*mm;
		//delete solidVest,logicVest,physVest,solidFilledVest;
		solidFilledVest = new G4Ellipsoid("VestFilledBox", fPhantomSize.x()/2 + 4.5*cm+thick,
						   fPhantomSize.y()/2 + 7.5*cm+thick,
						   fPhantomSize.z()/2 + 23.*cm+thick, 
						   -5.*cm,
						   fPhantomSize.z()/2-19.0*cm);
		
		G4cout << "ellipsoid dimensions " <<
		  fPhantomSize.x()/2 + 4.5*cm << " "
			<< fPhantomSize.y()/2 + 7.5*cm << " "
			<< fPhantomSize.z()/2 + 23.*cm << " " 
			<< -5.0*cm << " "
			<< fPhantomSize.z()/2-19.0*cm << G4endl;
				   
		solidVest = new G4SubtractionSolid("VestBox",solidFilledVest,containerSolid);						   
		logicVest = new G4LogicalVolume(solidVest,materials->GetMaterial("G4_WATER"),"logicVest");
		physVest = new G4PVPlacement(0, shift, logicVest, "physVest", logicAir, false, 0, true);
		//physVest = new G4PVPlacement(0, shift, logicVest, "physVest", logicWorld, false, 0, true);
		
		G4cout << "Volume of the vest = " << solidVest->GetCubicVolume() << G4endl;
	}

	G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void MyGeometry::SetCSVBodies(G4String csvBodies_)
{
	csvBodies = csvBodies_;
}


void MyGeometry::DefineCommands()
{
	fMessenger = new G4GenericMessenger(this, "/SIM/scoring/","Set scoring phantom");
	auto& phantomCmd = fMessenger->DeclareMethod("phantom",&MyGeometry::SetHumanPhantom,"Folder for results");
	phantomCmd.SetParameterName("phantom", true);
	phantomCmd.SetDefaultValue("IcruSphere");
	
	fMessenger = new G4GenericMessenger(this, "/SIM/scoring/","Set IHab module");
	auto& moduleCmd = fMessenger->DeclareMethod("putModule",&MyGeometry::SetModule,"Folder for results");
	moduleCmd.SetParameterName("putModule", true);
	moduleCmd.SetDefaultValue("0");
	
	
	fMessenger = new G4GenericMessenger(this, "/SIM/scoring/","Set CSV Bodies");
	auto& csvBodiesCmd = fMessenger->DeclareMethod("csvBodies",&MyGeometry::SetCSVBodies,"set csv bodies");
	csvBodiesCmd.SetParameterName("csvBodies", true);
	csvBodiesCmd.SetDefaultValue("../scene/bodies.csv");
	

	fMessenger = new G4GenericMessenger(this, "/SIM/geometry/vest/","Set thickness of vest");
	auto& thickVestCmd = fMessenger->DeclareMethod("thick",&MyGeometry::SetThickVest,"set vest thickness");
	thickVestCmd.SetParameterName("thick", true);
	thickVestCmd.SetDefaultValue("4");
}

