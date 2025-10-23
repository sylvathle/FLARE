#include "generator.hh"
#include "Randomize.hh"
#include "ions.hh"

#include "G4AnalysisManager.hh"

#include <math.h> 



MyPrimaryGenerator::MyPrimaryGenerator()
{
	fParticleGun = new G4ParticleGun(1);

	// Surface of beam in m2
	beamSurface = 4.0 * PI_ * rsource*rsource/m/m;

	iNbin = 100;
	ilogemin = 1;
	ilogemax = 5;

	ilogemin_gen = 1;
	ilogemax_gen = 5;

	primaryParticle = "proton";

	ions = getIons();
	
	ionTable = G4IonTable::GetIonTable();
	ionTable->CreateAllIon();
	/*G4cout << "Ions" << G4endl;
	for (G4int i=0;i<1000;i++)
	{
		G4cout << i << " " << ionTable->GetParticle(i)->GetParticleName() <<  G4endl;
	}
	G4cout << "ParticleTable" << G4endl;
	particleTable = G4ParticleTable::GetParticleTable();
	for (G4int i=0;i<510;i++)
	{
		G4cout << i << " " << particleTable->GetParticleName(i) << G4endl;
	}*/
	
	biasRndm = new G4SPSRandomGenerator();

	// Important scaling factor that rescale the angular distribution to 1/(4pi)
	// A factor pi must be introduced due to the cosine distribution over the surface of the 
	// sphere, then dividing by 4pi to normalize to isotropic distribution = 1/4
	solidAngleFactor = 0.25; 

	
	//for (int i=0;i<iNbin;i++) {Npart.push_back(0);}
	
	//sampleSize = 1e5;
	DefineCommands();

	std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
	auto duration = now.time_since_epoch();	
	auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);

	G4Random::setTheEngine(new CLHEP::RanecuEngine());
	G4long seed = nanoseconds.count()%10000000; //1;//time(NULL);
	G4Random::setTheSeed(seed);

	//nevents = 0;
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
	delete fMess;
}


void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
	//G4ParticleDefinition* particleDef = ionTable->GetIon(ions[primaryParticle].getZ(),ions[primaryParticle].getA(),0*keV);
	if (ions.find(primaryParticle) != ions.end()) 
	{
		G4cout << "Is an ion" << G4endl;
		particleDef = ionTable->GetIon(ions[primaryParticle].getZ(),ions[primaryParticle].getA(),0*keV);
	}
	else {particleDef = particleTable->FindParticle(primaryParticle);}
	
	
	G4double thetapos = CLHEP::RandFlat::shoot(-PI_,PI_);

	G4double posz0 = CLHEP::RandFlat::shoot(-rsource,rsource);
	G4double r = sqrt(rsource*rsource-posz0*posz0);
	G4double posx0 = r * cos(thetapos);
	G4double posy0 = r * sin(thetapos);

	
	G4double pow_energy = CLHEP::RandFlat::shoot(ilogemin_gen,ilogemax_gen);
	G4double energy = particleDef->GetBaryonNumber()*pow(10,pow_energy)*MeV;

	G4ThreeVector pos(posx0,posy0,posz0);

	G4ThreeVector mom = GenMomentum(pos);
	G4int iebin = iNbin*(pow_energy-ilogemin)/(ilogemax-ilogemin);	
	if (Npart.find(iebin) == Npart.end()) {Npart[iebin] = 1;} 
	else {Npart[iebin]++;}

	G4cout << particleDef->GetParticleName() << " " << iebin << " " << Npart[iebin] << G4endl;
	//nevents++;

	fParticleGun->SetParticlePosition(pos);
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleEnergy(energy);
	fParticleGun->SetParticleDefinition(particleDef);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}


G4ThreeVector MyPrimaryGenerator::GenMomentum(G4ThreeVector pos)
{
	// Normalize pos
	G4double normpos = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	G4double normposx = pos[0]/normpos;
	G4double normposy = pos[1]/normpos;
	G4double normposz = pos[2]/normpos;

	//Calculate vparallel
	G4double sizevpar = -sqrt(1-CLHEP::RandFlat::shoot(0.0,1.0));
	G4double xpar = normposx * sizevpar;
	G4double ypar = normposy * sizevpar;
	G4double zpar = normposz * sizevpar;

	
	// Calculate vperpendicular
	// First we need a base in the perpendicular plane
	// First base
	G4double normb1 = sqrt(zpar*zpar+xpar*xpar);
	G4double b1x = zpar/normb1;
	G4double b1y = 0.0;
	G4double b1z = -xpar/normb1;

	//Second base
	G4double b2x = normposy*b1z - normposz*b1y;
	G4double b2y = normposz*b1x - normposx*b1z;
	G4double b2z = normposx*b1y - normposy*b1x;

	
	G4double sizevperp = sqrt(1.0 - sizevpar*sizevpar);	


	G4double theta = CLHEP::RandFlat::shoot(-PI_,PI_);
	G4double mx0rot = sizevperp * (b1x * cos(theta) + b2x * sin(theta)) + xpar;
	G4double my0rot = sizevperp * (b1y * cos(theta) + b2y * sin(theta)) + ypar;
	G4double mz0rot = sizevperp * (b1z * cos(theta) + b2z * sin(theta)) + zpar;

	return G4ThreeVector(mx0rot,my0rot,mz0rot);

}

	
void MyPrimaryGenerator::SetPrimary(G4String particle_)
{
	primaryParticle = particle_;
}

void MyPrimaryGenerator::SetMinkE(G4double minlogkE)
{
	ilogemin_gen = minlogkE;
}

void MyPrimaryGenerator::SetMaxkE(G4double maxlogkE)
{
	ilogemax_gen = maxlogkE;
}


void MyPrimaryGenerator::SetBeamRadius(G4double r)
{
	rsource = r;
	beamSurface = 4.0 * PI_ * rsource*rsource/m/m;
}

void MyPrimaryGenerator::DefineCommands()
{
	//fMess = new G4GenericMessenger(this, "/SIM/scoring/","Set sample size");
	//auto& samplesizeCmd = fMess->DeclareMethod("sampleSize",&MyPrimaryGenerator::SetSampleSize,"Set size of sample");
	//samplesizeCmd.SetParameterName("sampleSize", true);
	//samplesizeCmd.SetDefaultValue("1e5");
	
	fMess = new G4GenericMessenger(this, "/SIM/scoring/","Set radius beam");
	auto& radiusbeamCmd = fMess->DeclareMethodWithUnit("radbeam","mm",&MyPrimaryGenerator::SetBeamRadius,"Set radius beam");
	radiusbeamCmd.SetParameterName("radbeam", true);
	radiusbeamCmd.SetDefaultValue("3004*mm");

	fMess = new G4GenericMessenger(this, "/SIM/generate/","Set primary");
	auto& partRatioCmd = fMess->DeclareMethod("particle",&MyPrimaryGenerator::SetPrimary,"Set primary");
	partRatioCmd.SetParameterName("particle", true);
	partRatioCmd.SetDefaultValue("H");

	fMess = new G4GenericMessenger(this, "/SIM/generate/","Set sample size");
	auto& logeminCmd = fMess->DeclareMethod("logminkE",&MyPrimaryGenerator::SetMinkE,"Set log of min kE of primary (per nucleon)");
	logeminCmd.SetParameterName("logminkE", true);
	logeminCmd.SetDefaultValue("1");

	fMess = new G4GenericMessenger(this, "/SIM/generate/","Set sample size");
	auto& logemaxCmd = fMess->DeclareMethod("logmaxkE",&MyPrimaryGenerator::SetMaxkE,"Set log of max kE of primary (per nucleon)");
	logemaxCmd.SetParameterName("logmaxkE", true);
	logeminCmd.SetDefaultValue("5");
}







