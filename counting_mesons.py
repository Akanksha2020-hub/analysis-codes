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
# Radius of jets (0.4, 1.0, 1.4) :
R = 1.4
# Events selection (pT in GeV)
pT_min_jet1 = 500
pT_min_jet2 = 500
eta_max = 2.5
M_PI = 3.14
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

#defining variables


# Get pointer to branches used in this analysis
# R-jet branches : 04, 08, 10, 12, 14
R_jet = str(int(R*10))
if R<1.0:
        R_jet = '0' + R_jet

# Getting the required branches from Delphes ROOT file.
branchJet = treeReader.UseBranch("ParticleFlowJet%s"%R_jet)
branchMET = treeReader.UseBranch("MissingET")
branchtrack = treeReader.UseBranch("Track")
branchParticle = treeReader.UseBranch("GenParticle")

hist1pion=ROOT.TH1F("diagonal_pions", "Number of diagonal pions in the Particle list", 100, 0.0, 100)
hist1rho=ROOT.TH1F("diagonal_rhos", "Number of diagonal rhos in the Particle list", 100, 0.0, 100)
hist2pion=ROOT.TH1F("offdiagonal_pions", "Number of off-diagonal pions in the Particle list", 100, 0.0, 100)
hist2rho=ROOT.TH1F("offdiagonal_rhos", "Number of off-diagonal rhos in the Particle list", 100, 0.0, 100)
hist3pion=ROOT.TH1F("total_pions", "Number of diagonal+off diagonal pions in the Particle list", 100, 0.0, 100)
hist3rho=ROOT.TH1F("total_rhos", "Number of off-diagonal rhos in the Particle list", 100, 0.0, 100)

for entry in range(0, numberOfEntries):
    treeReader.ReadEntry(entry)
    particlearray = branchParticle.GetEntries()
    nop1 = 0
    nor1 = 0
    nop2 = 0
    nor2 = 0
    
    for i in range (0,particlearray)
        pid = branchParticle.At(i).PID
    
        if abs(pid) == 4900111
            nop1 = nop1 + 1
        else:
            exit
        
        if abs(pid) == 4900113
            nor1 = nor1 + 1
        else:
            exit
        
        if abs(pid) == 4900211
            nop2 = nop2 + 1
        else:
            exit
        
        if abs(pid) == 4900213
            nor2 = nor2 + 1
        else:
            exit
        
    #Fill histograms
    hist1pion.Fill(nop1)
    hist1rho.Fill(nor1)
    hist2pion.Fill(nop2)
    hist2rho.Fill(nor2)
    hist3pion.Fill(nop1+nop2)
    hist3rho.Fill(nor1+nor2)
    counter=counter+1
        
#Printing number of accepted events           
print("counter=",counter)
integral=hist1pion.Integral(0,-1)
print("integral=",integral)
        
#Normalizing the histograms
if hist1pion.GetSumw2N()==0:
    hist1rho.Sumw2(True)
   
if hist1rho.GetSumw2N()==0:
    hist1rho.Sumw2(True)
        
if hist2pion.GetSumw2N()==0:
    hist2pion.Sumw2(True)
    
if hist2rho.GetSumw2N()==0:
    hist2rho.Sumw2(True)
        
if hist3pion.GetSumw2N()==0:
    hist3pion.Sumw2(True)
                       
if hist3rho.GetSumw2N()==0:
    hist3rho.Sumw2(True)
               
hist1pion.Scale(1./hist1pion.Integral()) 
hist1rho.Scale(1./hist1rho.Integral()) 
hist2pion.Scale(1./hist2pion.Integral()) 
hist2rho.Scale(1./hist2rho.Integral()) 
hist3pion.Scale(1./hist3pion.Integral())       
hist3rho.Scale(1./hist3rho.Integral())  

#Creating a list and saving the histograms to the list.

histlist = ROOT.TList()

histlist.Add(hist1pion)
histlist.Add(hist1rho)
histlist.Add(hist2pion)      
histlist.Add(hist2rho)
histlist.Add(hist3pion)
histlist.Add(hist3rho)

outputFile = inputFile[:-5] + "_counting_mesons_R14.root"
rootFile = ROOT.TFile(outputFile, "RECREATE")
histlist.Write()
rootFile.Close()


























