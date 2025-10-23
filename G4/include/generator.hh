#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

#include "G4ParticleGun.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleMomentum.hh"
#include "G4GenericMessenger.hh"
#include "G4SPSRandomGenerator.hh"

#include "G4RunManager.hh"
#include "geometry.hh"

#include "ions.hh"


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
		
		G4double GetGeometryFactor() const 
		{
			return solidAngleFactor*beamSurface;
		}
		G4String GetIonName() const { return primaryParticle;}

		// Return the flux for an ions Z with an arbitrary energy E
		std::map<G4int,G4int> GetNGenerated() const {return Npart;}
		G4int GetNGenerated(G4int i) {return Npart[i];}
		
		//G4int GetSampleSize() const {return sampleSize;}
		//G4int nevents;
		
		
	private :
	
		G4GenericMessenger *fMess;
		
		G4ParticleGun *fParticleGun;
  		G4double rsource, beamSurface; //, lowPosZ;

		G4double solidAngleFactor;

		G4SPSRandomGenerator *biasRndm;
		G4ParticleDefinition *particleDef;
		G4ParticleTable *particleTable;

		std::map<G4String,Ion> ions;
		G4IonTable *ionTable;
		G4String primaryParticle;

		//G4int sampleSize;

		G4int iNbin;
		G4double ilogemin,ilogemax;
		G4double ilogemin_gen,ilogemax_gen;

		std::map<G4int,G4int> Npart;

		G4ThreeVector GenMomentum(G4ThreeVector pos);
		void SetMinkE(G4double);
		void SetMaxkE(G4double);
		void SetPrimary(G4String);
		//void SetSampleSize(G4int s) {sampleSize=s;}
		void SetBeamRadius(G4double r);
		void DefineCommands();


	
};


#endif
