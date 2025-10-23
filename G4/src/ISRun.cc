//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// ISRun.cc
//
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

#include "ISRun.hh"


ISRun::ISRun()
:G4Run()
{
	//id_event=0;
	//id_entry=0;
	
	minlogE = 1.0;
	maxlogE = 5.0;
 	NbinsE = 100;

	const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	ISmass = detectorConstruction->GetLogicIS()->GetMass()/kg;
	
	generator = static_cast<const MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	mission_factor = generator->GetGeometryFactor()/joule;
	//n_event_record = generator->GetSampleSize();
}

ISRun::~ISRun()
{
	//delete generator;
 	fEDEMap.clear();
 	fDoseMap.clear();
	//Nparts.clear();
}

void ISRun::RecordEvent(const G4Event* event)
{
	// Hits collections
	//
	G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	G4int eventToBeProcessed = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
	
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) 
	{
	G4cout << "BYEEEE" << G4endl;	
	return;
	}

	//++id_event;

 	auto  EDE_id = G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/EDE");
  	
  	G4int energyBin = G4int((log10(event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy())-minlogE)/(maxlogE-minlogE)*NbinsE);

	auto* EDEMap = static_cast<G4THitsMap<HitInfo>*>(HCE->GetHC(EDE_id));
	// sum up the energy deposition and the square of it
	for (auto itr : *EDEMap->GetMap()) {
		fEDEMap[energyBin][0]  += itr.second->GetEDE(); //sum
		fDoseMap[energyBin][0]  += itr.second->GetDose();  //sum
	}

	//G4cout << n_event_record << G4endl;
	
	if (eventID+1==eventToBeProcessed)
	{
		RecordRunDoses();
		//id_event=0;
	}
}

void ISRun::RecordRunDoses()
{
	//G4cout << "RecordEvent" << G4endl;
	G4AnalysisManager *man = G4AnalysisManager::Instance();
	//G4cout << "RecordEvent" << G4endl;

	//G4int iter(0);//,id(0);
	G4double fac;
	G4double ede(0.0),HT(0.0),dose(0.0);
	G4int energyBin(0);
	//G4cout << "RecordEvent" << G4endl;

	//for (auto itr : massMap) {
	//G4String nameEntry = G4String(std::to_string(itr.first)+G4String("_")+G4String(particleName));
	for (auto imap : fEDEMap) {
		//G4cout << "RecordEvent" << G4endl;
		energyBin = imap.first;
		//G4cout << energyBin << " " << id_entry << G4endl;
			
		fac = mission_factor/ISmass*1e3; // #/m2/s/sr * m2 * s/y * sr / kg * J/# to mSv
		HT = fac*fEDEMap[energyBin][0];
		dose = fac*fDoseMap[energyBin][0];

		// Don't record the data if no energy deposited
		if ((HT==0) && (dose==0)){continue;}
		
		man->FillNtupleIColumn(0,0,id_entry);
		man->FillNtupleIColumn(0,1,energyBin);
		man->FillNtupleDColumn(0,2,HT);
		man->FillNtupleDColumn(0,3,dose);
			
		man->AddNtupleRow(0);
	}
	//++iter;
	//id_event=0;
	//++id_entry;

}

/*void ISRun::Merge(const G4Run* run)
{
	 // merge the data from each thread
	 EMAP localMap = static_cast<const ISRun*>(run)->fEDEMap;
	 G4cout << "Merging " << G4endl;
	 for(auto itr : localMap){
		 fEDEMap[itr.first]  += itr.second;
        	 //fEDEMap[itr.first].second += itr.second.second;
		}

	G4Run::Merge(run);
}*/




