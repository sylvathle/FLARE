#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

#include "G4ParticleGun.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleMomentum.hh"
#include "G4GenericMessenger.hh"
#include "G4SPSRandomGenerator.hh"

#include "G4RunManager.hh"
#include "geometry.hh"

#include "ions.hh"

//#include "functions.hh"

#include <vector>
#include <map>
#include <random>
#include <iostream>
#include <cmath>

#include <chrono>

//#include "ions.hh"
//#include "GCRFlux.hh"

#define PI_ 3.1415926535926535



class G4GeneralParticleSource;

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
	public :
		MyPrimaryGenerator();
		~MyPrimaryGenerator();
		
		virtual void GeneratePrimaries(G4Event*);
		G4double GetBeamSurface() const {return beamSurface;}
		G4double GetRSource() const {return rsource;}
		
		
		//G4double GetTotalFlux() const {return factorSphere*total_particles*beamSurface*yearToSecond*solidAngle;}
		//G4double GetTotalFlux() const {return listIntegratedFluxParticleKind[0]*beamSurface*yearToSecond*solidAngle;}
		//G4double GetTotalParticleNumber(G4int Z) const {return gcrFlux->GetParticleFlux(Z);}
		G4double GetMissionDuration() const {return yearToSecond;}
		G4double GetMissionFactor() const 
		{
			G4cout << "In GetMissionFactor " << beamSurface << G4endl;
			return factorSphere*beamSurface*solidAngle;
		}
		G4String GetIonName(G4int i) const { return listIons[i];}
		G4int GetNIons() const { return listIons.size();}
		//std::vector<G4String> GetParticleList() const {return ionsName;}
		//std::vector<G4double> GetGCRParticlesWeights() const {return gcrFlux->GetlistParticleWeight();}

		// Return the flux for an ions Z with an arbitrary energy E
		//G4double GetEnergyFlux(G4int Z, G4double E) const {return gcrFlux->GetEnergyFlux(Z,E);};
		//G4int GetNGenerated(G4int Z) const {return Npart[Z-1];}
		
		G4int GetSampleSize() const {return sampleSize;}
		
		
	private :
	
		G4GenericMessenger *fMess;
		
		G4ParticleGun *fParticleGun;
		G4bool halfSphere;
  		G4double rsource, beamSurface, lowPosZ,factorSphere;

		G4double yearToSecond,solidAngle;

		G4SPSRandomGenerator *biasRndm;
		G4ParticleTable *particleTable;

		//const GCRFlux *gcrFlux;
		//const G4String fluxFile;

		std::map<G4String,Ion> ions;
		G4IonTable *ionTable;

		G4int idParticle,nevent,nrejected,sampleSize;

		//std::vector<G4ParticleDefinition *> listParticles;
		//std::map<G4String,SBG4SPSEneDistribution *> listEneGenerator;

		std::vector<G4double> weightToUse;
		std::vector<G4double> ratioPart;
		std::vector<G4String> listIons;
		//std::vector<G4int> Npart;

		G4ThreeVector GenMomentum(G4ThreeVector pos);
		void SetParticleRatio(G4double,G4double);
		void SetSampleSize(G4int s) {sampleSize=s;}
		void SetBeamRadius(G4double r);
		void DefineCommands();


	
};

//G4double calculateTotalFlux(const std::string& filename);


#endif
