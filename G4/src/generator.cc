#include "generator.hh"
#include "Randomize.hh"
#include "ions.hh"

#include "G4AnalysisManager.hh"

#include <math.h> 



// Function to pick an index from the vector based on probabilities
int pickIndex(const std::vector<double>& probabilities) {
    // Calculate cumulative probabilities
    std::vector<double> cumulativeProbabilities(probabilities.size());
    double sum = 0.0;
    double factor = 1000000.0;
    for (size_t i = 0; i < probabilities.size(); ++i) {
        sum += probabilities[i]*factor;
        cumulativeProbabilities[i] = sum;
	//G4cout << i << " " << cumulativeProbabilities[i] << G4endl;
    }

    // Generate a random number between 0 and 1
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, factor);
    double randNum = dis(gen);

    // Find the index corresponding to the segment containing the random number
    for (size_t i = 0; i < cumulativeProbabilities.size(); ++i) {
	//G4cout << i << " " << randNum << " " << cumulativeProbabilities[i] << " " << cumulativeProbabilities.size() << G4endl;
        if (randNum < cumulativeProbabilities[i]) {
		//G4cout << "i = " << i << G4endl;
            return i;
        }
    }

    // In case of rounding errors, return the last index
    return probabilities.size() - 1;
}


//MyPrimaryGenerator::MyPrimaryGenerator():gcrFlux(new GCRFlux())
MyPrimaryGenerator::MyPrimaryGenerator()//:gcrFlux()
{
	fParticleGun = new G4ParticleGun(1);
	//fSources = new G4GeneralParticleSource();

	halfSphere = false;
	//halfSphere = true;

	nevent = 0;
	nrejected = 0;

	//rsource = 2.1*m;
	//rsource = 3.5*m;
	
	//rsource = static_cast<const MyGeometry*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction()->GetRSource());
	//const MyGeometry *detectorConstruction = static_cast<const MyGeometry*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	//rsource = detectorConstruction->GetRModule();
	
	//G4cout << "rsource = " << rsource << G4endl;
	
	beamSurface = 4.0 * PI_ * rsource*rsource;
	lowPosZ = -rsource;
	factorSphere = 1.0; //Case GCR is taken from all directions

	// If situation of the Lunar surface, the starting surface is half sphere, the mimal z position is 0, 
	//   and we are only considering half of the possible origins (the other being blocked by the Moon).
	if (halfSphere) 
	{
		beamSurface/=2.0;
		lowPosZ=0.0;
		factorSphere = 1.0;//0.5*.7213;
	}

	iNbin = 100; //bins.GetNbins();
	ilogemin = 1; //logMeV //bins.GetMinE(); //1 eV
	ilogemax = 5; //logMeV //bins.GetMaxE(); //100 GeV

	ilogemin_gen = 1; //logMeV //bins.GetMinE(); //1 eV
	ilogemax_gen = 5; //logMeV //bins.GetMaxE(); //100 GeV
	
	//particleTable = G4ParticleTable::GetParticleTable();
	ionTable = G4IonTable::GetIonTable();
	ionTable->CreateAllIon();

	//G4double excitEnergy = 0*keV;

	//G4ParticleDefinition* particleDef = ionTable->GetIon(26,56,excitEnergy);

	//G4cout << particleDef->GetParticleName() << G4endl;

	//ionTable->DumpTable();
	
	biasRndm = new G4SPSRandomGenerator();
	yearToSecond = 24.0 * 3600.0;
	solidAngle = PI_;

	//listIons.push_back("proton");
	//listIons.push_back("alpha");

	//G4int n_ions_kind = listIons.size();
	ions = getIons();
	
	//for (int i=0;i<n_ions_kind;i++) {Npart.push_back(0);} 
	//for (int i=0;i<n_ions_kind;i++) {ratioPart.push_back(G4double(1.0));} 
	//for (int i=0;i<n_ions_kind;i++) {weightToUse.push_back(G4double(1.0/n_ions_kind));} 
	//for (int i=0;i<n_ions_kind;i++) {SetParticleRatio(i+1,0);} 
	//SetParticleRatio(1,1);
	//for (int i=0;i<n_ions_kind;i++) {G4cout << weightToUse[i] << G4endl;}
	//weightToUse[0] = 1.0; 

	for (int i=0;i<iNbin;i++) {Npart.push_back(0);}
	
	sampleSize = 1e5;
	DefineCommands();


	//G4Random::setTheEngine(new CLHEP::RanecuEngine());
	//G4long seed = time(NULL);
	//G4Random::setTheSeed(seed);
	
	//G4cout << "Generator " << G4endl;
	
	std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
	auto duration = now.time_since_epoch();	
	auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);


	G4Random::setTheEngine(new CLHEP::RanecuEngine());
	G4long seed = nanoseconds.count()%10000000; //1;//time(NULL);
	G4Random::setTheSeed(seed);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
	delete fMess;
}


void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
	//idParticle = pickIndex(listParticleWeight);

	//idParticle = 0;
	//idParticle = pickIndex(weightToUse);
	//while (idParticle!=0 && idParticle!=1 && idParticle!=25) idParticle = pickIndex(listParticleWeight);
	//G4cout << idParticle << G4endl;
	//idParticle = 0;
	//G4String particleName = listIons[idParticle];
	
	
	//if (idParticle>7){G4cout << "Z = " << idParticle << " " << particleName << G4endl;}
	
	//G4ParticleDefinition* particleDef = ionTable->GetIon(ions[particleName].getZ(),ions[particleName].getA(),0*keV);
	//G4ParticleDefinition* particleDef = ionTable->GetIon(gcrFlux->GetZ(particleName),gcrFlux->GetA(particleName),0*keV);
	
	G4ParticleDefinition* particleDef = ionTable->GetIon(ions[primaryParticle].getZ(),ions[primaryParticle].getA(),0*keV);
	
	//G4double energy = listEneGenerator[particleName]->GenerateOne(particleDef);
	//energy = energy / 100.0;
	//G4cout << energy << G4endl;

	G4double thetapos = CLHEP::RandFlat::shoot(-PI_,PI_);

	G4double posz0 = CLHEP::RandFlat::shoot(-rsource,rsource);
	G4double r = sqrt(rsource*rsource-posz0*posz0);
	G4double posx0 = r * cos(thetapos);
	G4double posy0 = r * sin(thetapos);
	
	//posz0 = 2.2*m;
	//posy0 = -2.5*m;
	//posx0 = -.3*m;
	
	//energy = 0.5*GeV;
	//G4double pow_energy = CLHEP::RandFlat::shoot(1,2.6);
	G4double pow_energy = CLHEP::RandFlat::shoot(ilogemin_gen,ilogemax_gen);
	G4double energy = particleDef->GetBaryonNumber()*pow(10,pow_energy)*MeV;
	//energy = 10*MeV;	

	//G4cout << particleDef->GetParticleName() << " " << energy << G4endl;
	
	//G4cout << idParticle << " " << particleName << " " << energy << G4endl;
	//energy = 0.5*GeV;
	
	G4ThreeVector pos(posx0,posy0,posz0);

	//G4double momx = -posx0+CLHEP::RandFlat::shoot(-0.1,0.1)*m;
	//G4double momy = -posy0+CLHEP::RandFlat::shoot(-0.1,0.1)*m;
	//G4double momz = -posz0+CLHEP::RandFlat::shoot(-0.1,0.1)*m+0.3*m;

	G4ThreeVector mom = GenMomentum(pos);
	//G4double logprimkE = log10(energy);	
	G4int iebin = iNbin*(pow_energy-ilogemin)/(ilogemax-ilogemin);	
	Npart[iebin]++;
	//pos = G4ThreeVector(1.0,0.0,0.0);
	//mom = G4ThreeVector(momx,momy,momz);
	//G4cout << (G4double) (nrejected)/(G4double)(nevent) << G4endl;
	//G4cout << "GenParticles " << G4endl;

	//Npart[G4int(particleDef->GetPDGCharge())-1]++;

	//G4cout << particleName << " " << gcrFlux->GetZ(particleName) << " " << energy << G4endl;
	

	fParticleGun->SetParticlePosition(pos);
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleEnergy(energy);
	fParticleGun->SetParticleDefinition(particleDef);
	fParticleGun->GeneratePrimaryVertex(anEvent);
}



// Function generating a cosine distribution from the source point
G4ThreeVector MyPrimaryGenerator::GenMomentum(G4ThreeVector pos)
{
	// Normalize pos
	G4double normpos = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	G4double normposx = pos[0]/normpos;
	G4double normposy = pos[1]/normpos;
	G4double normposz = pos[2]/normpos;
	//G4cout << normposx << " " << normposy << " " << normposz << G4endl;

	G4double u = CLHEP::RandFlat::shoot(0.0,1.0);
	G4double v = CLHEP::RandFlat::shoot(0.0,1.0);

	//G4cout << "u = " << u << "  v=" << v << G4endl;

	// theta and phi in spherical coordinates
	G4double theta = asin(u);
	G4double phi = 2 * PI_ * v;

	//G4cout << theta << " " << phi << G4endl;

	// generated vector along the z axis
	G4double x = sin(theta)*cos(phi);
	G4double y = sin(theta)*sin(phi);
	G4double z = cos(theta);

	// Calculate vparallel
	G4double xpar = normposy;
	G4double ypar = -normposx;
	G4double zpar = 0.0;
	G4double normpar = sqrt(xpar*xpar + ypar*ypar + zpar*zpar);
	if (normpar==0)
	{
		xpar = -normposz;
		ypar = 0;
		zpar = normposx;
	}
	normpar = sqrt(xpar*xpar + ypar*ypar + zpar*zpar);
	xpar = xpar/normpar;
	ypar = ypar/normpar;
	zpar = zpar/normpar;

	// Calculate vperp
	G4double xper = normposy*zpar - normposz*ypar;
	G4double yper = normposz*xpar - normposx*zpar;
	G4double zper = normposx*ypar - normposy*xpar;
	
	G4double normper = sqrt(xper*xper + yper*yper + zper*zper);
	xper = xper/normper;
	yper = yper/normper;
	zper = zper/normper;

	//New vector
	G4double xrot = x * xpar + y * xper + z * normposx;
	G4double yrot = x * ypar + y * yper + z * normposy;
	G4double zrot = x * zpar + y * zper + z * normposz;
	G4double normrot = sqrt(xrot*xrot + yrot*yrot + zrot*zrot);
	xrot = -xrot/normrot;
	yrot = -yrot/normrot;
	zrot = -zrot/normrot;
	
	return G4ThreeVector(xrot,yrot,zrot);
}

	
void MyPrimaryGenerator::SetPrimary(G4String particle_)
{
//	ratioPart[G4int(Z)-1] = rat;
	primaryParticle = particle_;
	
	// Normalize
	//G4double tot = 0.;
	//for (G4int i=0;i<weightToUse.size();i++){tot+=ratioPart[i];}
	//if (tot==0) return;
	//for (G4int i=0;i<weightToUse.size();i++){weightToUse[i]=ratioPart[i]/tot;}
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
	//G4cout << "BeamSurafce " << beamSurface << G4endl;
	lowPosZ = -rsource;

	// If situation of the Lunar surface, the starting surface is half sphere, the mimal z position is 0, 
	//   and we are only considering half of the possible origins (the other being blocked by the Moon).
	if (halfSphere) 
	{
		beamSurface/=2.0;
		lowPosZ=0.0;
		factorSphere = 1.0; //0.5*.7213;
	}
	
	//G4cout << "SetBeamRadius  " << rsource << G4endl;
}

void MyPrimaryGenerator::DefineCommands()
{
	fMess = new G4GenericMessenger(this, "/SIM/scoring/","Set sample size");
	auto& samplesizeCmd = fMess->DeclareMethod("sampleSize",&MyPrimaryGenerator::SetSampleSize,"Set size of sample");
	samplesizeCmd.SetParameterName("sampleSize", true);
	samplesizeCmd.SetDefaultValue("1e5");
	
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







