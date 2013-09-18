#! /usr/bin/env python
import ROOT as r
from TableMaker import *

period    = 'Combined'
selection = 'muGamma'
dataSet   = 'Muons'
doZJets   = True

LUMIDATA = 0.
if selection in ['muon', 'muGamma']:
    if period in ['Combined', '2011A']:
        LUMIDATA += .2151 
        LUMIDATA += .9302 
        LUMIDATA += .3709 
        LUMIDATA += .6630 
        #LUMIDATA *= 1.0514
    if period in ['Combined', '2011B']:
        LUMIDATA += 2.511
elif selection in ['electron', 'eGamma']:
    if period in ['Combined', '2011A']:
        LUMIDATA += .2151 
        LUMIDATA += .7982 
        LUMIDATA += .3132 
        LUMIDATA += .6332 
        #LUMIDATA *= 1.053
    if period in ['Combined', '2011B']:
        LUMIDATA += 2.482

bTable = TableMaker('../HiggsAnalyzer/histos/higgsHistograms_'+dataSet+period+'.root', LUMIDATA)

samples = ['tW', 'ttbar', 'GluGluWW', 'WWJets', 'WZJets', 'ZZJets']
if selection in ['muGamma', 'eGamma']:
    samples.append('PhotonJets')
elif selection == 'muEG':
    samples.append('DYToTauTau')
elif doZJets:
    samples.append('ZJets')
elif not doZJets:
    if selection == 'muon':
        samples.append('DYToMuMu')
        samples.append('DYToTauTau')
    if selection == 'electron':
        samples.append('DYToEE')
        samples.append('DYToTauTau')


bTable.AddDatasets(samples)
bTable.AddDatasets('DATA')

bTable._rowList     = ['0 b-jet', '1 b-jet', '2 b-jet', '3 b-jet']
if selection in ['eGamma', 'muGamma']:
    bTable._columnList  = ['ZZJets', 'WZJets', 'WWJets', 'ttbar', 'tW', 'PhotonJets', 'BG', 'DATA'] 
if selection in ['electron', 'muon']:
    bTable._columnList  = ['ZZJets', 'WZJets', 'WWJets', 'ttbar', 'tW', 'ZJets', 'BG', 'DATA'] 

histDict = bTable.GetHistDict('Jet', 'bJetMultNoVeto')
bTable.PrintTable(histDict, True, True, startBin = 1)
