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
	id_event=0;
	id_entry=0;
	
	minlogE = -5.0;
	maxlogE = 5.0;
 	NbinsE = 10;

	const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	ISmass = detectorConstruction->GetLogicIS()->GetMass()/kg;
	//fTetData = detectorConstruction->GetICRP145phantom();
	
	//for (auto itr : massMap) {massMap[itr.first]=massMap[itr.first]/kg;}
	//organsGrouped = GetOrgansGroup("organsInfo.csv");
	//mapGroupedOrgans = generateNumberedMap(organsGrouped);
	//N_groups = GetNGroups(mapGroupedOrgans);
	//wT = GetwT("organsInfo.csv");
  //	csvFile = std::ofstream("results.csv");
	
	//const MyPrimaryGenerator *generator = static_cast<const MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	generator = static_cast<const MyPrimaryGenerator*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	//mission_factor = generator->GetTotalParticleKindNumber();
	//mission_factor = generator->GetMissionFactor()/joule;
	n_event_record = generator->GetSampleSize();
	
	//G4cout << "n_event_record " << n_event_record << G4endl;
	//rsource = generator->GetRSource()/m;
	//GCRParticleWeights = generator->GetGCRParticlesWeights();
	//total_flux = generator->GetTotalFlux()/joule;
	//for (G4int i=0;i<generator->GetNIons();i++)
	//{
	//	//G4cout << i << " " << generator->GetTotalParticleNumber(i) << G4endl;
	//	GCRParticleWeights.push_back(mission_factor*generator->GetTotalParticleNumber(i));
	//}
	//G4cout << "factor = " << mission_factor << G4endl;
	for (G4int Z=1;Z<=28;Z++) {Nparts.push_back(0);}
	
}

ISRun::~ISRun()
{
	//delete generator;
 	fEDEMap.clear();
 	fDoseMap.clear();
	Nparts.clear();
}

void ISRun::RecordEvent(const G4Event* event)
{
	// Hits collections
	//
	
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) 
	{
	G4cout << "BYEEEE" << G4endl;	
	return;
	}

	++id_event;

 	auto  EDE_id = G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/EDE");
  	
  	G4int Z = (G4int)event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetPDGCharge();
	//G4double E = event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
	//G4double d = generator->GetEnergyFlux(Z,E);
	//d = 1.0;
	
	//G4cout << Z << " " << E << " " << d << G4endl;

	//G4cout << "ISR " << rsource << " " << z << " " << (G4double) (totalAngle) / (G4double)(n_events) << G4endl;

  	G4int energyBin = G4int((log10(event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy())-minlogE)/(maxlogE-minlogE)*NbinsE);

	//G4cout << "Z = " << Z << G4endl;
	
	//Z=0;
	//Nparts[Z-1]++;
	Nparts[Z-1]++;
  	//G4cout << "Z = " << Z << G4endl;

	//G4cout << id_event << " " << Z << G4endl;
  	
	auto* EDEMap = static_cast<G4THitsMap<HitInfo>*>(HCE->GetHC(EDE_id));
	// sum up the energy deposition and the square of it
	for (auto itr : *EDEMap->GetMap()) {
		//G4String nameEntry = G4String(std::to_string(itr.first))+G4String("_")+G4String(std::to_string(Z)));
		fEDEMap[energyBin][Z]  += itr.second->GetEDE(); //sum
		fDoseMap[energyBin][Z]  += itr.second->GetDose();  //sum
		//G4cout << id_event << " Z=" << Z << " " << itr.first << " " << fEDEMap[Z] << " " << fDoseMap[Z] << " " << Nparts[Z-1] << " " << itr.second->GetEDE() << G4endl;
		
		//if (fEDEMap[id_entry]!=0.0) {G4cout << "fedemap " << id_entry << " " << fEDEMap[id_entry] << G4endl;}
		//G4cout << id_entry << " " << id_event << " " << itr.second->GetEDE() << " " << fEDEMap[0] << G4endl;
	}

	
	//G4cout << "SSSSSSSSSS  " << fEDEMap.size() << G4endl;
	if (id_event==n_event_record)
	{
		G4AnalysisManager *man = G4AnalysisManager::Instance();

		//G4int totpart=0;
		
		for (G4int iZ=0;iZ<=27;iZ++){
			man->FillNtupleIColumn(0,iZ,Nparts[iZ]);
			//totpart+=Nparts[iZ];
		}
		//G4cout << id_event << " " << id_entry << " " << "totalparts = " << totpart << G4endl;
		man->AddNtupleRow(0);
		G4int iter(1);//,id(0);
		//G4int i_organ = 0;
		G4double fac;
		G4double ede(0.0),HT(0.0),dose(0.0);

		//for (auto itr : massMap) {
		//G4String nameEntry = G4String(std::to_string(itr.first)+G4String("_")+G4String(particleName));
		man->FillNtupleIColumn(iter,0,id_entry);
		man->FillNtupleIColumn(iter,1,0);
		for (G4int iZ=0;iZ<=27;iZ++)
		{

			//G4int iZ = Z;
			if (Nparts[iZ]>0)
			{
				fac = 1.0/ISmass/Nparts[iZ]*1e3; // #/m2/s/sr * m2 * s/y * sr / kg * J/# to mSv
				//G4String Zstring = G4String(std::to_string(iZ+1));
				HT = fac*fEDEMap[energyBin][iZ];
				dose = fac*fDoseMap[energyBin][iZ];
				//if (iZ==1) { G4cout << id_entry << " " << HT << " " << Nparts[iZ] << G4endl;}
			}
			else {HT=0.0; dose=0.0;}
			//G4cout << iter << " " << itr.first << " " << itr.second << " " << id_entry << " " << N_organ << " " << HT << G4endl;
			man->FillNtupleDColumn(iter,2*(iZ+1),HT);
			man->FillNtupleDColumn(iter,2*(iZ+1)+1,dose);
			
			//G4cout << itr.first << " " << wT[itr.first] << G4endl;
			if (iZ>=0) {ede += HT;}
			//tot_dose += HT/N_organ;
		}
		man->AddNtupleRow(1);
		//++i_organ;
		//}
		++iter;
		man->FillNtupleDColumn(iter,0,ede);
		//man->FillNtupleDColumn(iter,1,tot_dose*mission_factor*1e3);
		man->AddNtupleRow(iter);	
		
		//fEDEMap.clear();
		//fDoseMap.clear();
		for (G4int iZ=0;iZ<28;iZ++) 
		{
			Nparts[iZ]=0;
			//G4String Zstring = G4String(std::to_string(iZ+1));
			fEDEMap[energyBin][iZ+1]=0;
			fDoseMap[energyBin][iZ+1]=0;
		}

	
		id_event=0;
		++id_entry;
	}
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




