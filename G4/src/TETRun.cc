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
// TETRun.cc
//
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

#include "TETRun.hh"


TETRun::TETRun()
:G4Run()
{
	id_event=0;
	id_entry=0;
	
	minlogE = 1.0;
	maxlogE = 5.0;
 	NbinsE = 100;
	
	const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	fTetData = detectorConstruction->GetICRP145phantom();
	massMap = fTetData->GetMassMap();
	
	N_organ= massMap.size();
	
	for (auto itr : massMap) 
	{
		massMap[itr.first]=massMap[itr.first]/kg;
	}
	const MyPrimaryGenerator *generator = static_cast<const MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	//G4cout << std::setprecision(20) << 1.0/joule << G4endl;
	//G4cout << std::setprecision(20) << 1.0*joule << G4endl;
	mission_factor = generator->GetGeometryFactor()/joule;
	n_event_record = generator->GetSampleSize();
}

TETRun::~TETRun()
{
 	fEDEMap.clear();
 	fDoseMap.clear();
	Nparts.clear();
}

void TETRun::RecordEvent(const G4Event* event)
{
	++id_event;
	
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

 	auto  EDE_id = G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/EDE");
  	
  	G4int Z = (G4int)event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetPDGCharge();
  	G4double A = (G4double)event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetBaryonNumber();
	G4double energy = event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy()/A;
  	G4int energyBin = G4int((log10(energy)-minlogE)/(maxlogE-minlogE)*NbinsE);
	G4double energyBinMin = pow(10,energyBin/(G4double)NbinsE*(maxlogE-minlogE)+minlogE);
	G4double energyBinMax = pow(10,(energyBin+1)/(G4double)NbinsE*(maxlogE-minlogE)+minlogE);
	G4double ebinsize = energyBinMax - energyBinMin;

	Z = 1;
	Nparts[energyBin][Z-1]++;
  	
	auto* EDEMap = static_cast<G4THitsMap<HitInfo>*>(HCE->GetHC(EDE_id));
	// sum up the energy deposition and the square of it
	for (auto itr : *EDEMap->GetMap()) {
		G4int id = itr.first+Z*100000; 
		
		fEDEMap[energyBin][id]  += itr.second->GetEDE(); //sum
		fDoseMap[energyBin][id]  += itr.second->GetDose();  //sum
	}
	
	if (id_event==n_event_record)
	{
		G4AnalysisManager *man = G4AnalysisManager::Instance();
		
		G4int iter(0), id(0);
		G4double fac;
		G4double HT(0.0),dose(0.0);

		for (auto imap : fEDEMap) {
			energyBin = imap.first;
			for (auto itr : massMap) {
				man->FillNtupleIColumn(iter,0,itr.first);
				man->FillNtupleIColumn(iter,1,energyBin);
				//mission_factor = (Joules).(m2) / massMap (kg)  ->>> mGy.m2
				fac = mission_factor/massMap[itr.first]*1e3;
				for (G4int iZ=0;iZ<=1;iZ++)
				{
					id = itr.first+(iZ+1)*100000;
					if (Nparts[energyBin][iZ]<=0)
					{
						continue;
					}
					
					HT = fac*fEDEMap[energyBin][id];
					dose = fac*fDoseMap[energyBin][id];
					man->FillNtupleDColumn(iter,2*iZ+2,HT);
					man->FillNtupleDColumn(iter,2*iZ+3,dose);
				}
				man->AddNtupleRow(0);
			}
		}
		++iter;
		
		fEDEMap.clear();
		fDoseMap.clear();
		Nparts.clear();
	
		id_event=0;
		++id_entry;
	}
}

/*void TETRun::Merge(const G4Run* run)
{
	 // merge the data from each thread
	 EMAP localMap = static_cast<const TETRun*>(run)->fEDEMap;
	 G4cout << "Merging " << G4endl;
	 for(auto itr : localMap){
		 fEDEMap[itr.first]  += itr.second;
        	 //fEDEMap[itr.first].second += itr.second.second;
		}

	G4Run::Merge(run);
}*/




