#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4CSGSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4NistManager.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4SystemOfUnits.hh"
#include "G4MultiFunctionalDetector.hh"

//#include "SBG4PSFlatSurfaceFlux.hh"
#include "SBG4PSSphereSurfaceFlux.hh"

#include "G4VPrimitiveScorer.hh"
#include "G4RunMessenger.hh"
#include "G4GenericMessenger.hh"
//#include "G4UIcmdWithADoubleAndUnit.hh"

//#include "detector.hh"
//#include "SBMultiFunctionalDetector.hh"
//#include "shieldmaterial.hh"

#include <sstream>


#include "G4SDManager.hh"



#include "TETModelImport.hh"
#include "TETParameterisation.hh"
#include "TETPSEnergyDeposit.hh"
#include "ISPSEnergyDeposit.hh"

#include "materials.hh"
#include "generator.hh"

#define USE_CADMESH_TETGEN TRUE


class MyGeometry : public G4VUserDetectorConstruction
{
	public:		
		MyGeometry();
		MyGeometry(G4UIExecutive *ui_);
		~MyGeometry();
		
		G4LogicalVolume *GetLogicIS() const { return logicPhantom; }

		virtual G4VPhysicalVolume *Construct();
		TETModelImport * GetICRP145phantom() const {return fTetData;}
		G4String GetPhantomType() const {return phantomType;}
		//G4double GetRModule() const {return radiusModule+thickModule;}
	
	private:
		G4UIExecutive *ui;
		
		// rPhantom stands for the radius of the Human Phantom
		//void DefineMaterials();
		void SetHumanPhantom(G4String);
		void SetModule(G4double);

		void SetThickVest(G4double);
		void SetCSVBodies(G4String);
		

		//void SetRegolithDome(G4double inner_radius, G4double outer_radius);
		void DefineCommands();
		
		//G4Material* GetMaterial(G4String mat);
		//G4Material* GetRegoAtDepth(G4double);
		
		
	protected:
		//G4Material *worldMat, *icruSphereMat, *RegoBrick, *AluMat, *EAC1A, *PLA, *LiqMethane, *LiqHydrogene, *RegoAp17, *Cr2O3, *MnO, *P2O5;

		Materials *materials;
	
		G4double rPhantom ;
		//G4LogicalVolume *fScoringVolume;
		G4GenericMessenger *fMessenger;

		G4double heightModule;
		G4double radiusModule;
		G4double thickModule;		
		
		G4ThreeVector shift;
		G4Sphere *solidPhantom;
		//G4Ellipsoid* containerSolid;
		G4VSolid* containerSolid;
		G4Box *solidWorld;
		G4Sphere *solidAir;
		TETModelImport *fTetData;
		G4LogicalVolume *logicWorld, *logicAir,*logicPhantom, *fTetLogic, *fContainer_logic;
		G4VPhysicalVolume *physWorld, *physAir, *physPhantom;
		
		// Vest related objects
		G4Ellipsoid* solidFilledVest;
		G4VSolid *solidVest;
		G4LogicalVolume *logicVest;
		G4VPhysicalVolume *physVest;

		// Module related objetcs
		//G4Tubs *solidModule, *solidModuleUp, *solidModuleDown;
		G4Sphere *solidModule;
		G4LogicalVolume *logicModule, *logicModuleUp, *logicModuleDown;
		G4VPhysicalVolume *physModule,  *physModuleUp, *physModuleDown;

		G4double thickVest;
		
                G4String phantomType;
                G4String csvBodies;

                // rPhantom stands for the radius of the Human Phantom
                //G4double zWorld;
                G4ThreeVector      fPhantomSize;
                G4ThreeVector      fPhantomBoxMin, fPhantomBoxMax;
	        G4int              fNOfTetrahedrons;

                virtual void ConstructSDandField();

};



#endif 
