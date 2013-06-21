#! /usr/bin/env python
import sys, os, pickle
import ROOT as r
from TableMaker import *

dataPath  = '~/nobackup/HZZ/archive/2011/Dec06/higgsHistograms_'
#dataPath  = '../HiggsAnalyzer/histos/higgsHistograms_'
period    = 'Combined'
selection = 'electron'
dataSet   = 'Electrons'
doZJets   = True  
doPresel  = True  
doSignal  = False

LUMIDATA = 0.
if selection == 'muon' or selection == 'muGamma':
    if period in ['2011A', 'Combined']:
        LUMIDATA += .2151 
        LUMIDATA += .9302 
        LUMIDATA += .3709 
        LUMIDATA += .6630 
    if period in ['2011B', 'Combined']:
        LUMIDATA += 2.511
elif selection == 'electron' or selection == 'eGamma':
    if period in ['2011A', 'Combined']:
        LUMIDATA += .2151 
        LUMIDATA += .7982 
        LUMIDATA += .3132 
        LUMIDATA += .6332 
    if period in ['2011B', 'Combined']:
        LUMIDATA += 2.482


samples = ['tW', 'ttbar', 'GluGluWW', 'WWJets', 'WZJets', 'ZZJets']

if selection == 'muGamma' or selection == 'eGamma':
    samples.append('PhotonJets')
elif not doZJets:
    if selection == 'muon':
        samples.append('DYToMuMu')
        samples.append('DYToTauTau')
    if selection == 'electron':
        samples.append('DYToEE')
        samples.append('DYToTauTau')
else:
    samples.append('ZJets')

yieldTable = TableMaker(dataPath+dataSet+period+'.root', LUMIDATA, 'electron', delimiter = '&')

  ## PhotonJets workaround ##

yieldTable._scaleDict['PhotonJets'] = 0.5

  ## PRESELECTION ##

if doPresel:
    yieldTable._rowList     = ['single lepton', 'HLT selection/weights', 'Good PV', 'Data quality', 'Cosmics rejection', 'Dilepton', 'Z mass', 'Z p<sub>T</sub>', '3rd lepton veto', 'b-jet veto']
    if selection in ['muGamma', 'eGamma']:
        yieldTable._columnList  = ['ZZJets', 'WZJets', 'WWJets', 'ttbar', 'tW', 'PhotonJets', 'BG', 'DATA'] 
    else:
        yieldTable._columnList  = ['ZZJets', 'WZJets', 'WWJets', 'ttbar', 'tW', 'ZJets', 'BG', 'DATA'] 

    yieldTable.AddDatasets(samples)
    yieldTable.AddDatasets('DATA')

    histDict = yieldTable.GetHistDict('Misc', 'acceptanceByCut')
    yieldTable.PrintTable(histDict, True, False, startBin = 1)

  ## SIGNAL ##
if doSignal:
    yieldTable._rowList      = ['250', '300', '350', '400', '450', '500', '550', '600']
    if selection in ['muGamma', 'eGamma']:
        #yieldTable._columnList  = ['ZZJets', 'WZJets', 'WWJets', 'ttbar', 'tW', 'PhotonJets', 'BG', 'Signal', 'DATA'] 
        yieldTable._columnList  = ['ZZJets', 'WZJets', 'WWJets', 'top', 'PhotonJets', 'BG', 'Signal', 'DATA'] 
    else:
        yieldTable._columnList  = ['ZZJets', 'WZJets', 'WWJets', 'ttbar', 'tW', 'ZJets', 'BG', 'Signal', 'DATA'] 

    yieldTable.AddDatasets(samples, Clear=True)
    yieldTable.AddDatasets(['HZZ250', 'HZZ300', 'HZZ350', 'HZZ400', 'HZZ450', 'HZZ500', 'HZZ550', 'HZZ600'])
    yieldTable.AddDatasets(['HWW250', 'HWW300', 'HWW350', 'HWW400', 'HWW450', 'HWW500', 'HWW550', 'HWW600'])
    yieldTable.AddDatasets(['VBF250', 'VBF300', 'VBF350', 'VBF400', 'VBF450', 'VBF500', 'VBF550', 'VBF600'])
    yieldTable.AddDatasets('DATA')

    histDict = yieldTable.GetHistDict2D('Misc', 'nEventsByHMass', 4) 
    yieldTable.PrintTable(histDict, True, True, startBin = 1)
