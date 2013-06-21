#! /usr/bin/env python
from FitMaster import *

### For PU weights ###
fitter = FitMaster(inputFile = '../data/PU/Cert_177718_178078_7TeV_PromptReco_Collisons11_JSON.pileup_v2.root', inHist = 'pileup', dataName = 'pileUp_run2011B')
fitter.DoFit('[0]*(TMath::Poisson(x, [1]))', 0., 35., [10000, 5])

### For photon weights ###
#fitter = FitMaster(inputFile = '../data/Reweight.root', inHist = 'h1_eGammaPtWeight2011B', dataName = 'ptWeight_electron')
#fitter.DoFit('[0]+[1]*x^[2]', 170., 420., [-2., 1.8, 0.05])

#fitter = FitMaster(inputFile = '../data/higgsHistograms_Muons2011A.root', inHist = 'h1_diLeptonQt_DATA', dataName = 'dileptonQt',  directory = 'Lepton')
#fitter.DoFit('([0]/[1])*(x**2/[1])**([2]-1)*exp(-1*(x**2-[3])/(2*[1]))', 55., 200., [100000., 1.43, 2.5, 1.])
