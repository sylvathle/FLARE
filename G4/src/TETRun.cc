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
		//G4cout << itr.first << " " << massMap[itr.first]/kg << G4endl;
	}
	//organsGrouped = GetOrgansGroup("organsInfo.csv");
	//mapGroupedOrgans = generateNumberedMap(organsGrouped);
	//N_groups = GetNGroups(mapGroupedOrgans);
	//wT = GetwT("../scripts/organsInfo.csv");
  //	csvFile = std::ofstream("results.csv");
	
	const MyPrimaryGenerator *generator = static_cast<const MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	//mission_factor = generator->GetTotalParticleKindNumber();
	mission_factor = generator->GetMissionFactor()/joule;
	n_event_record = generator->GetSampleSize();
	//GCRParticleWeights = generator->GetGCRParticlesWeights();
	//for (G4int i=0;i<GCRParticleWeights.size();i++) {GCRParticleWeights[i]=generator->GetTotalParticleNumber(i)/joule;}
	//for (G4int i=0;i<generator->GetNIons();i++)
	//{
	//	//G4cout << i << " " << generator->GetTotalParticleNumber(i) << G4endl;
	//	GCRParticleWeights.push_back(mission_factor*generator->GetTotalParticleNumber(i));
	//}
	//G4cout << "factor = " << mission_factor << G4endl;
	//for (G4int Z=1;Z<=1;Z++) {Nparts.push_back(0);}

}

TETRun::~TETRun()
{
 	fEDEMap.clear();
 	fDoseMap.clear();
	Nparts.clear();
}

void TETRun::RecordEvent(const G4Event* event)
{
	// Hits collections
	//
	++id_event;
	
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

 	auto  EDE_id = G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/EDE");
  	
  	G4int Z = (G4int)event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetPDGCharge();
  	G4double A = (G4double)event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetBaryonNumber();
  	//G4cout << A << " " << event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy()/A <<  G4endl;
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
		//G4String nameEntry = G4String(std::to_string(itr.first))+G4String("_")+G4String(std::to_string(Z)));
		G4int id = itr.first+Z*100000; 
		//G4cout << energyBin << " " << Z << " " << id << G4endl;
		
		fEDEMap[energyBin][id]  += itr.second->GetEDE(); //sum
		fDoseMap[energyBin][id]  += itr.second->GetDose();  //sum
		//G4cout << energyBin << " " << id << " " << fEDEMap[energyBin][id] << G4endl;
		//if (fEDEMap[id_entry]!=0.0) {G4cout << "fedemap " << id_entry << " " << fEDEMap[id_entry] << G4endl;}
	}
	
	if (id_event==n_event_record)
	{
		G4AnalysisManager *man = G4AnalysisManager::Instance();
		//for (auto ip : Nparts)
		//{
			//man->FillNtupleIColumn(0,iZ,Nparts[iZ]);
			//man->FillNtupleIColumn(0,0,ip.first);
			//for (auto iz : ip.second)
			//{
				//man->FillNtupleIColumn(0,iz.first,ip.first)
				//G4cout << ip.first << " " << iz.first << " " << iz.second << G4endl;
			//	man->FillNtupleIColumn(0,iz.first+1,iz.second);
			//}
			//man->AddNtupleRow(0);
		//}
		//for (G4int iZ=0;iZ<=0;iZ++){man->FillNtupleIColumn(0,iZ,Nparts[iZ]);}
		
		G4int iter(0), id(0);
		//G4int i_organ = 0;
		G4double fac;
		G4double HT(0.0),dose(0.0);

		for (auto imap : fEDEMap) {
			energyBin = imap.first;
			energyBinMin = pow(10,energyBin/(G4double)NbinsE*(maxlogE-minlogE)+minlogE);
			energyBinMax = pow(10,(energyBin+1)/(G4double)NbinsE*(maxlogE-minlogE)+minlogE);
			ebinsize = energyBinMax - energyBinMin;
			for (auto itr : massMap) {
				//G4String nameEntry = G4String(std::to_string(itr.first)+G4String("_")+G4String(particleName));
				//G4cout << id_entry << " " << itr.first << G4endl;
				//man->FillNtupleIColumn(iter,0,id_entry);
				//fac = mission_factor/massMap[itr.first]*1e3;
				man->FillNtupleIColumn(iter,0,itr.first);
				man->FillNtupleIColumn(iter,1,energyBin);
				//G4cout << "mission_factor " << mission_factor << G4endl;
				//G4cout << itr.first << " Mass " << massMap[itr.first] << G4endl;
				fac = mission_factor/massMap[itr.first]*1e3*ebinsize; //mission_factor (m2.sr) / massMap (kg) / N  ->>> mGy.m2.sr/MeV / particle
				for (G4int iZ=0;iZ<=1;iZ++)
				{
					id = itr.first+(iZ+1)*100000;
					if (Nparts[energyBin][iZ]<=0)
					{
						//G4cout << fEDEMap[energyBin][id] << " " << energyBin << " " << iZ << " " << Nparts[energyBin][iZ] << " " << GCRParticleWeights[iZ] << G4endl;
						//fac = GCRParticleWeights[iZ]/massMap[itr.first]/Nparts[energyBin][iZ]*1e3;
						//HT = fEDEMap[energyBin][id]/massMap[itr.first]/Nparts[energyBin][iZ]*1e3/joule;
						//dose = fDoseMap[energyBin][id]/massMap[itr.first]/Nparts[energyBin][iZ]*1e3/joule;
						continue;
					}
					
					//G4cout << iZ << " " << energyBinMin << " " << energyBinMax << " " << fDoseMap[energyBin][id] << "  " << fEDEMap[energyBin][id] << "  " << itr.first << "  " << massMap[itr.first] <<  "  " << ebinsize << G4endl;
					HT = fac*fEDEMap[energyBin][id];
					dose = fac*fDoseMap[energyBin][id];
					//G4cout << iter << " " << itr.first << " " << itr.second << " " << id_entry << " " << N_organ << " " << HT << G4endl;
					//man->FillNtupleIColumn(iter,3*iZ+2,Nparts[energyBin][iZ]);
					man->FillNtupleDColumn(iter,2*iZ+2,HT);
					man->FillNtupleDColumn(iter,2*iZ+3,dose);
					//G4cout << iter << " " << 2*(iZ+1) << " " << Nparts[energyBin][iZ] << " " << HT << " " << dose << G4endl;
					//G4cout << itr.first << " " << wT[itr.first] << G4endl;
					//ede += HT*wT[itr.first];
					//tot_dose += HT/N_organ;
				}
				man->AddNtupleRow(0);
				//++i_organ;
			}
		}
		++iter;
		//man->FillNtupleDColumn(iter,0,ede);
		//man->FillNtupleDColumn(iter,1,tot_dose*mission_factor*1e3);
		//man->AddNtupleRow(iter);	
		
		fEDEMap.clear();
		fDoseMap.clear();
		Nparts.clear();
		//for (G4int iZ=0;iZ<=0;iZ++) {Nparts[iZ]=0;}
	
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




