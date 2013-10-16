#! /usr/bin/env python
from array import array
from PlotProducer import *
import subprocess
import ROOT as r


def ratio_1D(ratioDict, path, inFile, outFile):

    for key,value in ratioDict.iteritems():

        print 'Producing ratio plot {0}'.format(key)

        h1_Numer = inFile.GetDirectory(path).Get('h1_{0}'.format(value[0]))
        h1_Denom = inFile.GetDirectory(path).Get('h1_{0}'.format(value[1]))

        g_Eff = r.TGraphAsymmErrors()
        g_Eff.Divide(h1_Numer, h1_Denom)
        g_Eff.SetName('g_{0}'.format(key))
        g_Eff.SetTitle('{0};p_{{T}};#varepsilon'.format(key))

        outFile.Add(g_Eff)


def ratio_2D(ratioDict, path, inFile, outFile):

    h2_Eff = []

    for key,value in ratioDict.iteritems():

        print 'Producing ratio plot {0}'.format(key)

        h2_Numer = inFile.GetDirectory(path).Get('h2_{0}'.format(value[0]))
        h2_Denom = inFile.GetDirectory(path).Get('h2_{0}'.format(value[1]))

        h2_Eff.append(r.TH2D('h2_{0}'.format(key), '{0};;'.format(key),
                         h2_Numer.GetNbinsX(), h2_Numer.GetXaxis().GetXmin(), h2_Numer.GetXaxis().GetXmax(),
                         h2_Numer.GetNbinsY(), h2_Numer.GetYaxis().GetXmin(), h2_Numer.GetYaxis().GetXmax()))
        h2_Eff[len(h2_Eff)-1].Divide(h2_Numer, h2_Denom)

        outFile.Add(h2_Eff[len(h2_Eff)-1])


if __name__ == '__main__':

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)


    ### For electron charge misID efficiencies ###
    #inFile  = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20131003_171426.root', 'OPEN')
    #outFile = r.TFile('data/electronQMisID.root', 'RECREATE')

    #eMisQDict = {
    #    'LeadElectronMisQ':('LeadElecQMisIDNumer', 'LeadElecQMisIDDenom'),
    #    'TrailingElectronMisQ':('TrailingElecQMisIDNumer', 'TrailingElecQMisIDDenom'),
    #    'DielectronMisQ':('DileptonQMisIDNumer', 'DileptonQMisIDDenom')}

    #ratio_2D(eMisQDict, 'inclusive/DATA_ELECTRON', inFile)

    ### For lepton fake rate measurement ###
    inFile  = r.TFile('fakeEstimator/fakeHistograms.root', 'OPEN')
    outFile = r.TFile('data/fakeRates.root', 'RECREATE')

    fakeDict1D = {
        'MuonFakePt_Even':('MuPassLepPt', 'MuProbeLepPt'),
        'MuonFakeEta_Even':('MuPassLepEta', 'MuProbeLepEta'),
        'MuonFakePt':('MuNumerPt', 'MuDenomPt'),
        'MuonFakeEta':('MuNumerEta', 'MuDenomEta'),
    }
    ratio_1D(fakeDict1D, 'inclusive/TEST', inFile, outFile)

    fakeDict2D = {
        'MuonFake2D':('MuNumer', 'MuDenom')#,
        #'ElectronFake2D':('EleNumer', 'EleDenom')
    }
    ratio_2D(fakeDict2D, 'inclusive/TEST', inFile, outFile)

    print 'Done!!!'
    outFile.Write()
    outFile.Close()

