#include "TotalFlux.hh"

TotalFlux::TotalFlux()
{
	// Initialize data collection
	//bins = Bins();
	Nbin = 3000; //bins.GetNbins();
	logemin = -9; //logMeV //bins.GetMinE(); //100 KeV
	logemax = 6; //logMeV //bins.GetMaxE(); //100 GeV

	nIncident = 0;

	//for (G4int i=0;i<Nbin;i++) {nIncident.push_back(0);}
	
	//for (G4int i=0;i<=Nbin;i++)
	//{
	//	edgesLog.push_back( logemin + (logemax-logemin)*i/Nbin );
	//	edges.push_back(pow(10,logemin + (logemax-logemin)*i/Nbin));
	//}	
}

TotalFlux::TotalFlux(G4String csvSource)
{
	this->FirstUse();
}


TotalFlux::~TotalFlux()
{}

void TotalFlux::FirstUse()
{
	// Initialize data collection
	Nbin = 3000;// bins.GetNbins();
	logemin = -9; //bins.GetMinE(); //100 KeV
	logemax = 6; //bins.GetMaxE(); //100 GeV
	
	//G4cout << "TotalFlux FIRST USE " << G4endl; 

	nIncident = 0;
	//for (G4int i=0;i<Nbin;i++) {nIncident.push_back(0);}

}

/*void TotalFlux::ToCSV(G4String fileName) const
{
	G4cout << "TOCSV 1" << G4endl;
	std::ofstream csvFile(fileName);
	if (csvFile.is_open())
	{
		csvFile << "particle,N_primary,iebin,iabin";
		for (int i=0;i<Nbin*Nangles;i++)
		{
			csvFile << ",obi_" << i ;
		}
		csvFile << "\n";



		std::map<G4int, ParticleSpectra>::iterator it;
		for (auto const& x : ofluxes)
		{
			ParticleSpectra spe = x.second; // string's value 
			for (int iebin=0;iebin<Nbin;iebin++)
			{
			
				//G4cout << "iabin " << iabin << G4endl;
				for (int iabin=0;iabin<Nangles;iabin++)
				{
					//G4cout << " iebin " << iebin << " iabin " << iabin << " index " << iebin*Nangles+iabin << G4endl;
					csvFile << x.first << "," << nIncident[iabin*Nbin+iebin] << "," << iebin << "," << iabin ;
					for (int oebin=0;oebin<Nbin;oebin++)
					{						
						for (int oabin=0;oabin<Nangles;oabin++)
						{
							//G4cout << iabin << " " << iebin << " " << oabin << " " << oebin << G4endl;
							 csvFile << ","  << spe.GetBin(iabin,iebin,oabin,oebin);
						}
					}
					csvFile << "\n";
				}
			}			
			
		}
	}
	csvFile.close();
}*/


void TotalFlux::IterIncident()
{
	//G4double logprimkE = log10(segkE);
	//if (logprimkE>=logemin)
	//{
	//	G4int oebin = Nbin*(logprimkE-logemin)/(logemax-logemin);
		//G4cout << "iterincident" << iebin<< " " << iabin << " " << nIncident[iabin*Nbin+iebin]  << G4endl;
		nIncident++;
		//G4cout << " nIncident " << nIncident << G4endl;
	//}
}


// Update total flux considering energy of primary, secondary and secondary identity
// If energies out of spectra boundaries, abort update and leave function
bool TotalFlux::Update(G4double secondarykE, G4String secondaryName)
{
	//std::cout << "Beggining Update function" << std::endl;
	//std::cout << secondarykE << " " << secondaryName << std::endl;

	//if (secondaryName=="nu_mu" || secondaryName=="nu_e" || secondaryName=="anti_nu_mu" || secondaryName=="anti_nu_e") {return false;}
	if (secondaryName!="neutron") {return false;}
	//if (secondaryName!="proton" && secondaryName!="neutron" && secondaryName!="gamma" && secondaryName!="pi+" && secondaryName!="pi-" && secondaryName!="e+" && secondaryName!="e-") {return false;}
	G4double logsegkE = log10(secondarykE);
	if (logsegkE<logemin) {return false;}
	
	//G4int iebin = Nbin*(logprimkE-logemin)/(logemax-logemin);
	//if (iebin<0) return false;
	//if (iebin>=Nbin)
	//{
	///	G4cout << "Warning ibin too high energy particle for the energy range provided, ignoring" << G4endl; 
	//	return false;
	//}


	G4int oebin = Nbin*(logsegkE-logemin)/(logemax-logemin);
	if (oebin<0) return false;
	if (oebin>=Nbin)
	{
		G4cout << "Warning ibin too high energy particle for the energy range provided, ignoring" << G4endl; 
		return false;
	}

	//std::cout << "obini defined " << obin << std::endl;
	//G4int iabin = G4int(primaryAng*Nangles);
	//if (primaryAng == 1.0) {iabin = Nangles-1;}
	//G4int oabin = G4int(secondaryAngle*Nangles);
	//if (secondaryAngle == 1.0) {oabin = Nangles-1;}
	
	G4String particleName(secondaryName);
	//G4cout << "TotalFlux " << particleName << " updateflussssss" << iebin << " " << oabin << G4endl;

	if (ofluxes.count(particleName)==0) {ofluxes[particleName]=ParticleSpectra(Nbin);}
	ofluxes[particleName].Update(oebin);
	
	//G4cout << "TotalFlux " <<  "FLUX UPDATED" << G4endl;
	
	//nIncident[ibin]++;

	//ofluxes[secondaryName].print();
	
	return true;
}

void TotalFlux::Print()
{

	std::map<G4int, ParticleSpectra>::iterator it;
	for (auto const& x : ofluxes)
	{
		std::cout << x.first << std::endl; 
		ParticleSpectra spe = x.second; // string's value 
		spe.print();
	}
	G4cout << "nIncident " << nIncident << G4endl;
}

