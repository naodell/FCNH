#! /usr/bin/env python
import subprocess, shlex, time, math, sys, os, pickle
from array import array

import ROOT as r
from GridSearch import *

#paramFile   = open('scripts/fcncParams.pkl', 'rb')
#scales      = pickle.load(paramFile)
#styles      = pickle.load(paramFile)
#combos      = pickle.load(paramFile)
#removes     = pickle.load(paramFile)
#categories  = pickle.load(paramFile)
#systematics = pickle.load(paramFile)

def add_systematic(hist, category, dataset):
    '''
    Add a systematic uncertainty to a histogram based the data sample and the
    category.
    '''

    if dataset not in systematics[category]:
        return hist

    for bin in range(hist.GetNbinsX()):
        entries = hist.GetBinContent(bin+1)
        if entries == 0:
            continue

        errorSq = pow(hist.GetBinError(bin+1), 2)
        for syst in systematics[category][dataset]:
            errorSq += pow(syst*entries, 2)

        error = math.sqrt(errorSq)
        #print self._category, dataset, entries, error, hist.GetBinError(bin+1)
        hist.SetBinError(bin+1, error)

    return hist


if __name__ == '__main__':

    if len(sys.argv) > 0:
        batch   = sys.argv[1]
    else:
        print 'A batch and ratio type must be specified.  Otherwise, do some hacking so this thing knows about your inputs.'
        exit()

    r.gStyle.SetOptStat(0)
    #categories  = ['ss_ee', 'ss_emu', 'ss_mumu', 'ss_inclusive']
    categories  = ['ss_inclusive']
    backgrounds = ['irr', 'Fakes', 'QFlips']
    datasets    = ['data_obs', 'fcnh'] + backgrounds

    dataDict = {}
    dataDict['data_obs']    = ['DATA_ELECTRON', 'DATA_MUON', 'DATA_MUEG'] # Observed
    dataDict['fcnh']        = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_ZZ_t', 'FCNC_ZZ_tbar', 'FCNC_TauTau_t', 'FCNC_TauTau_tbar'] # signal

    dataDict['irr']         = ['WZJets3LNu', 'ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau', 'ttZ', 'ttW', 'ttG'] # Irreducible backgrounds
    dataDict['Fakes']       = ['muFakes', 'eFakes', 'llFakes'] # Fakes
    dataDict['QFlips']      = ['QFlips'] # electron charge misID

    histFile   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_{0}.root'.format(batch), 'OPEN')
    #histFile   = r.TFile('fcncAnalysis/histos/fcncHistograms_cut1.root'.format(batch), 'OPEN')

    lumi        = 19.7e3
    variable    = 'BDT'

    # prepare output directory
    filePath = 'data/shapes'
    if not os.path.exists(filePath):
        os.system('mkdir -p '+filePath)
    elif len(os.listdir(filePath)) is not 0:
        os.system('rm -r {0}/*'.format(filePath))


    yields  = {}
    sumHists = []
    for category in categories:
        outfile = r.TFile('data/shapes/{0}_{1}_{2}.root'.format(category, variable, batch), 'RECREATE')
        yields[category] = CatData(category, datasets)

        # Get histograms for different samples 
        for dataset in datasets:
            sumHist = None 
            for i,sample in enumerate(dataDict[dataset]):
                histDir = histFile.GetDirectory(category + '/' + sample)
                if not histDir: continue

                hist = histDir.Get('h1_' + variable)
                if not hist: continue

                ### Do scaling of MC samples ###
                if dataset not in ['data_obs', 'QFlips', 'Fakes']:
                    yieldHist   = histFile.GetDirectory('inclusive/' + sample).Get('h1_YieldByCut')
                    nInit       = yieldHist.GetBinContent(1)
                    hist.Scale(scales['2012'][sample]*lumi/nInit)

                if dataset in ['Fakes', 'QFlips']:
                    hist = add_systematic(hist, category, sample)
                else:
                    hist = add_systematic(hist, category, 'Irreducible')

                if sumHist == None:
                    sumHist = hist.Clone()
                    sumHist.SetName(dataset)
                else:
                    sumHist.Add(hist)


            if sumHist:
                sumHists.append(sumHist)

                if dataset != 'data_obs': 
                    ### Make systematic shapes ###
                    histUp      = sumHist.Clone()
                    histDown    = sumHist.Clone()
                    histUp.SetName('{0}_Up'.format(dataset))
                    histDown.SetName('{0}_Down'.format(dataset))

                    for bin in range(histUp.GetNbinsX()):
                        histUp.SetBinContent(bin + 1, histUp.GetBinContent(bin+1)+histUp.GetBinError(bin+1))
                        histDown.SetBinContent(bin + 1, histDown.GetBinContent(bin+1)-histDown.GetBinError(bin+1))

                    sumHists.append(histUp)
                    sumHists.append(histDown)

                numEvents = sumHist.Integral()
                #print category,dataset,numEvents

                yields[category].add_data(dataset, numEvents)

        outfile.Write()
        outfile.Close()

    dataCard = open('{0}/{1}_{2}.txt'.format(filePath, variable, 'shape'), 'w')
    card_producer(yields, categories, backgrounds, dataCard)

