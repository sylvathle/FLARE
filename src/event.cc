#include "event.hh"
#include <csignal>

MyEventAction::MyEventAction(MyRunAction* myrun)
{
	myRun = myrun;
}

MyEventAction::~MyEventAction()
{}


void MyEventAction::BeginOfEventAction(const G4Event*)
{
	totalSurfFluxID = 1;
}



void MyEventAction::EndOfEventAction(const G4Event* evt)
{
	//G4String primaryName = evt->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetParticleName();
	G4double ikE = evt->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
	myRun->IterIncident(ikE);
	G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
	G4THitsMap<fluxInfo>* eventTotalSurfFlux = (G4THitsMap<fluxInfo>*)(HCE->GetHC(totalSurfFluxID));
	std::map<G4int, fluxInfo*>::iterator it;
	for (it = eventTotalSurfFlux->GetMap()->begin(); it != eventTotalSurfFlux->GetMap()->end(); it++)
	{
		myRun->UpdateFlux(ikE,it->second->kE,it->second->particleName);
		//G4cout << "event.cc  UpdateFLUX done" << G4endl;
	}
}
