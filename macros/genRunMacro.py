import sys

""" Scenario can be:
	0 - test
	1 - B2G-naked_freespace
	2 - B2G-naked
	3 - B2G-vest
	4 - ICRP-vest
	5 - ICRP-naked
"""

#os.environ["PHANTOM_PATH"] = "/geant4/v11.1.1/geant4-v11.1.1/examples/advanced/ICRP145_HumanPhantoms/build/ICRP145data/"

scenario = sys.argv[1]
#scenario = "ICRP-naked"

if "ICRP" in scenario: phantom = "ICRP145"
else: phantom = "BDRTOG4"


nsim = int(sys.argv[2])
#nsim = randrange(100000)

thicknessModule = int(sys.argv[3])
#thicknessModule = 4

particle = sys.argv[4]
logemin = sys.argv[5]
logemax = sys.argv[6]

innerRadiusModule = 1200
distSourceModule = 200

radiusSource = thicknessModule + innerRadiusModule + distSourceModule

fmacros_name = "../macros/run_"+scenario+"_"+particle+"_"+str(thicknessModule)+".mac"

fmacros = open(fmacros_name,'w')
fmacros.write("/run/verbose 0\n")

if phantom=="BDRTOG4":
  if "vest" in scenario:
    fmacros.write("/SIM/scoring/csvBodies ../scene/bodiesAndVest.csv\n")
  else:
    fmacros.write("/SIM/scoring/csvBodies ../scene/bodiesOnly.csv\n")

fmacros.write("/SIM/scoring/phantom "+phantom+"\n\n")

fmacros.write("/run/initialize\n\n")

fmacros.write("/SIM/scoring/putModule "+str(thicknessModule)+"\n")
fmacros.write("/SIM/scoring/sampleSize "+str(nsim)+"\n")
fmacros.write("/SIM/scoring/resDir "+"../results/"+scenario+"_"+str(thicknessModule)+"_"+particle+"\n")
#fmacros.write("/SIM/scoring/resDir "+"/data/results/"+particle+"/"+scenario+"_"+str(thicknessModule)+"\n")
fmacros.write("/SIM/scoring/radbeam " + str(radiusSource) + " mm\n\n")

fmacros.write("/SIM/generate/particle " + particle + "\n")
fmacros.write("/SIM/generate/logminkE " + logemin + "\n")
fmacros.write("/SIM/generate/logmaxkE " + logemax + "\n\n")

fmacros.write("/run/beamOn "+str(nsim)+"\n")

print (scenario,nsim,thicknessModule)
