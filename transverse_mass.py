command: python /path_of_code/transverse_mass.py /path_of_rootfile/name_of_rootfile.root

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

# Parameters
############################################
# Radius of jets (0.4, 0.8, 1.0, 1.2, 1.4) :
R = 1.4
# Events selection (pT in GeV)
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

# Book histograms
Nbins = 150
hmt = ROOT.TH1F("jet_met_mt" , "Transverse mass of jet + MET" , Nbins, 0.0 , 3000.0)

vec1 = ROOT.TLorentzVector()
vec2 = ROOT.TLorentzVector()
vec3 = ROOT.TLorentzVector()
vec4 = ROOT.TLorentzVector()
vec5 = ROOT.TLorentzVector()
vec6 = ROOT.TLorentzVector()
mt = 0
counter=0
integral=0
vecjet = ROOT.TLorentzVector()
vec = ROOT.TLorentzVector()
# Loop over all events
for entry in range(0, numberOfEntries):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)

    # If event contains at least 2 jet
    if branchJet.GetEntries() > 1:
        # Take the two leading jets
        jet1 = branchJet.At(0)
        jet2 = branchJet.At(1)
        
        if jet1.PT > pT_min_jet1 and np.abs(jet1.Eta) < eta_max and jet2.PT > pT_min_jet2 and np.abs(jet2.Eta) < eta_max :
            # Defining the TLorentz vector for leading jet
            vec1 = ROOT.TLorentzVector()
            vec1.SetPtEtaPhiM(jet1.PT , jet1.Eta , jet1.Phi , jet1.Mass )

            vec2 = ROOT.TLorentzVector()
            vec2.SetPtEtaPhiM(jet2.PT , jet2.Eta , jet2.Phi , jet2.Mass)

            px1 = jet1.PT*np.cos(jet1.Phi)
            py1 = jet1.PT*np.sin(jet1.Phi)
            pz1 = jet1.PT*np.sinh(jet1.Eta)
            energy1 = np.sqrt((jet1.Mass * jet1.Mass) + (jet1.PT*np.cosh(jet1.Eta) * jet1.PT*np.cosh(jet1.Eta)))
            vec3 = ROOT.TLorentzVector()
            vec3.SetPxPyPzE(px1 , py1 , pz1 , energy1)

            px2 = jet2.PT*np.cos(jet2.Phi)
            py2 = jet2.PT*np.sin(jet2.Phi)
            pz2 = jet2.PT*np.sinh(jet2.Eta)
            energy2 = np.sqrt((jet2.Mass * jet2.Mass) + (jet2.PT*np.cosh(jet2.Eta) * jet2.PT*np.cosh(jet2.Eta)))
            vec4 = ROOT.TLorentzVector()
            vec4.SetPxPyPzE(px2, py2 , pz2 , energy2)
            vecjet = ROOT.TLorentzVector()
            vecjet = vec3 + vec4
            counter=counter+1
            # Defining TLorentz vector for missing energy branch
            met1 = branchMET.At(0)
            m1 = met1.MET
            METPhi1 = met1.Phi
            METx1 = m1* np.cos(METPhi1)
            METy1 = m1* np.sin(METPhi1)

            vec5 = ROOT.TLorentzVector()
            vec5.SetPxPyPzE(METx1 , METy1 , 0 , m1)

            #Adding the vectors
            vec = ROOT.TLorentzVector()
            vec = vecjet + vec5
            mt = (vec).Mt()

            # Filling transverse mass histogram
            hmt.Fill(mt)

print("counter=",counter)
integral=hmt.Integral(0,-1)
print("integral=",integral)

#Normalizing the histogram
if hmt.GetSumw2N()==0:
    hmt.Sumw2(True)

hmt.Scale(1./hmt.Integral())
histlist = ROOT.TList()

histlist.Add(hmt)
outputFile = inputFile[:-5] + "_mt.root"
rootFile = ROOT.TFile(outputFile, "RECREATE")
histlist.Write()
rootFile.Close()


















