#! /usr/bin/env python
import subprocess, shlex, datetime, copy, math, pickle
from multiprocessing import Process
import ROOT as r
import PlotProducer as pp

paramFile   = open('scripts/fcncParams.pkl', 'rb')
scales      = pickle.load(paramFile)
styles      = pickle.load(paramFile)
combos      = pickle.load(paramFile)


if __name__ == '__main__':

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)

    fakeFile    = r.TFile('fakeEstimator/histos/20141130_234120.root', 'OPEN')

    datasets    = ['WJets', 'ZJets', 'ttbar', 'QCD']#, 'DATA']
    sampleSets  = {
                'DATA':['DATA_MUON', 'DATA_ELECTRON', 'DATA_MUEG'],
                'WJets':['W1JetsToLNu', 'W2JetsToLNu', 'W3JetsToLNu', 'W4JetsToLNu'],
                'ZJets':['ZJets_M-50', 'ZJets_M-10To50'],
                'ttbar':['ttbarHad', 'ttabarLep'],
                'QCD':['QCD_20_MU', 'QCD_20-30_EM', 'QCD_30-80_EM', 'QCD_80-170_EM', 'QCD_170-250_EM', 'QCD_250-350_EM', 'QCD_350_EM']
                }   

    canvas = r.TCanvas('canvas', 'canvas', 600, 450)

    legend = r.TLegend(0.65, 0.6, 0.89, 0.89)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    lumi        = 1.
    selection   = 'MC_Truth'
    variables   = ['MatchedMuonJetFlavor_inclusive']
    hists       = {}
    for dataset in datasets:
        hists[dataset] = 0.
        if dataset in sampleSets.keys():
            sampleHists = {}
            for sample in sampleSets[dataset]: 
                hYields = fakeFile.Get('inclusive/{0}/h1_YieldByCut'.format(sample))
                if not hYields: continue
                norm = 1.
                if dataset != 'DATA':
                    nInit   = hYields.GetBinContent(1)
                    norm    = lumi*scales['2012'][sample]/nInit
                
                for variable in variables:
                    hVar = fakeFile.Get('{0}/{1}/h1_{2}'.format(selection, sample, variable))
                    if not hVar: continue
                    hVar.Scale(1./norm)
                    hists[dataset] = hVar
                    
        else:
            hYields = fakeFile.Get('inclusive/{0}/h1_YieldByCut'.format(dataset))
            if not hYields: continue
            if dataset[:4] != 'DATA':
                nInit   = hYields.GetBinContent(1)
                norms[dataset] = lumi*scales['2012'][sample]/nInit

            for variable in variables:
                hVar = fakeFile.Get('{0}/{1}/h1_{2}'.format(selection, sample, variable))
                if not hVar: continue
                hVar.Scale(1./norm)
                hists[dataset] = hVar
                    

    print hists 

    #for var in variables:
    #    for dataset in datasets:
