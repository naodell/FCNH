#! /usr/bin/env python
import subprocess, shlex, time, math, sys, os, pickle
from array import array
import ROOT as r


if __name__ == '__main__':

    if len(sys.argv) > 0:
        batch   = sys.argv[1]
    else:
        print 'A batch and ratio type must be specified.  Otherwise, do some hacking so this thing knows about your inputs.'
        exit()

    r.gStyle.SetOptStat(0)
    categories  = ['ss_ee', 'ss_emu', 'ss_mumu', 'ss_inclusive']
    backgrounds = ['irr', 'Fakes', 'QFlips']
    datasets    = ['obs', 'signal'] + backgrounds

    dataDict = {}
    dataDict['obs']     = ['DATA_ELECTRON', 'DATA_MUON', 'DATA_MUEG'] # Observed
    dataDict['signal']  = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_ZZ_t', 'FCNC_ZZ_tbar', 'FCNC_TauTau_t', 'FCNC_TauTau_tbar'] # signal

    dataDict['irr']     = ['WZJets3LNu', 'ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau', 'ttZ', 'ttW', 'ttG'] # Irreducible backgrounds
    dataDict['Fakes']   = ['muFakes', 'eFakes', 'llFakes'] # Fakes
    dataDict['QFlips']  = ['QFlips'] # electron charge misID

    histFile   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_{0}.root'.format(batch), 'OPEN')
    #histFile   = r.TFile('fcncAnalysis/histos/fcncHistograms_cut1.root'.format(batch), 'OPEN')

    lumi        = 19.7
    paramFile   = open('scripts/fcncParams.pkl', 'rb')
    scales      = pickle.load(paramFile)
    variable    = 'BDT'

    sumHists = []
    for category in categories:
        outfile = r.TFile('data/shapes/{0}_{1}_{2}.root'.format(category, variable, batch), 'RECREATE')
        # Get histograms for different samples 
        for dataset in datasets:
            sumHist = None 
            for i,sample in enumerate(dataDict[dataset]):
                histDir = histFile.GetDirectory(category + '/' + sample)
                if not histDir: continue

                hist = histDir.Get('h1_' + variable)
                if not hist: continue

                ### Do scaling of MC samples ###
                if dataset not in ['obs', 'QFlips', 'Fakes']:
                    yieldHist   = histFile.GetDirectory('inclusive/' + sample).Get('h1_YieldByCut')
                    nInit       = yieldHist.GetBinContent(1)
                    hist.Scale(scales['2012'][sample]*lumi/nInit)

                if sumHist == None:
                    sumHist = hist.Clone()
                    sumHist.SetName(dataset)
                else:
                    sumHist.Add(hist)


            if sumHist:
                sumHists.append(sumHist)

                if dataset != 'obs': 
                    ### Make systematic shapes ###
                    histUp      = r.TH1D()
                    histDown    = r.TH1D()
                    histUp      = sumHist.Clone()
                    histDown    = sumHist.Clone()
                    histUp.SetName('{0}_Up'.format(dataset))
                    histDown.SetName('{0}_Down'.format(dataset))

                    for bin in range(histUp.GetNbinsX()):
                        histUp.SetBinContent(bin + 1, histUp.GetBinContent(bin)+histUp.GetBinError(bin))
                        histDown.SetBinContent(bin + 1, histDown.GetBinContent(bin)+histDown.GetBinError(bin))

                    sumHists.append(histUp)
                    sumHists.append(histDown)


        outfile.Write()
        outfile.Close()
