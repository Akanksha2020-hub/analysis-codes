#command: python /path_of_code/trackpt.py /path_of_rootfile/name_of_rootfile.root

import sys
import numpy as np
import ROOT
from array import array

try:
  input = raw_input
except:
  pass

if len(sys.argv) < 2:
  print(" Usage: Analysis_code/Jet_analysis.py /path/delphes_file.root")
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
        ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
        ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
        pass

############################################
# Radius of jets (0.4, 0.8, 1.0, 1.2, 1.4) :
R = 1.4
M_PI = 3.14
pT_min_jet1 = 500
pT_min_jet2 = 500
eta_max = 2.5
############################################

inputFile = sys.argv[1]
print("Input file :")
print(inputFile)

# Create chain of root trees

chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointer to branches used in this analysis
# R-jet branches : 04, 08, 10, 12, 14
R_jet = str(int(R*10))
if R<1.0:
        R_jet = '0' + R_jet

branchJet = treeReader.UseBranch("ParticleFlowJet%s"%R_jet)
branchMET = treeReader.UseBranch("MissingET")
branchtrack = treeReader.UseBranch("Track")

eta1 = 0
eta2 = 0
delEta = 0
counter=0
dEta = 0
dPhi = 0
averagept = 0
DeltaR = 0

Nbins = 20
histtrackpt = ROOT.TH2F("track_pt", "track PT vs Ntrk1", 20,0.0,60.0 , 20,0.0,160.0)
histdelEta = ROOT.TH1F("delta_eta", "delta_eta_jet1_and_jet2", 50,0.0,10.0)
# Loop over all events
for entry in range(0, numberOfEntries):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)

    # If event contains at least 2 jet
    if branchJet.GetEntries() > 1:
        # Take the two leading jets
        jet1 = branchJet.At(0)
        jet2 = branchJet.At(1)
        
        #if jet1.PT > pT_min_jet1 and np.abs(jet1.Eta) < eta_max and jet2.PT > pT_min_jet2 and np.abs(jet2.Eta) < eta_max :
        #print("ntrk1=",jet1.NCharged)
        eta1 = jet1.Eta
        eta2 = jet2.Eta
        delEta = jet1.Eta - jet2.Eta
        histdelEta.Fill(delEta)
        track = []
        selectedtrack = []
        trackpt = []

        for i in range(0 , branchtrack.GetEntries()):
             track.append(branchtrack.At(i))

        for j in range(0 , len(track)):
            dPhi = track[j].Phi - jet1.Phi
            if (dPhi  >  M_PI): dPhi -= 2*M_PI
            if (dPhi <= - M_PI): dPhi += 2*M_PI
            dEta = track[j].Eta - jet1.Eta
            DeltaR = np.sqrt(dEta*dEta + dPhi*dPhi)
            if DeltaR < 1.4 :
                selectedtrack.append(track[j])
        for k in range(0 , len(selectedtrack)):
            trackpt.append(selectedtrack[k].PT)

        averagept = np.average(trackpt)
        histtrackpt.Fill(averagept,jet1.NCharged)
            counter=counter+1
#Normalizing the histograms    
histdelEta.Scale(1./histdelEta.Integral())
histlist = ROOT.TList()
histlist.Add(histtrackpt)
histlist.Add(histdelEta)

outputFile = inputFile[:-5] + "_trackpt.root"
rootFile = ROOT.TFile(outputFile, "RECREATE")
histlist.Write()
rootFile.Close()


