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


def ratio_2D(ratioDict, path, inFile):

    outFile = r.TFile('electronQMisID.root', 'RECREATE')

    h2_Eff = []

    for key,value in ratioDict.iteritems():

        print 'Producing ratio plot {0}'.format(key)

        h2_Numer = inFile.GetDirectory(path).Get('h2_{0}'.format(value[0]))
        h2_Denom = inFile.GetDirectory(path).Get('h2_{0}'.format(value[1]))

        h2_Eff.append(r.TH2D('h2_{0}'.format(key), '{0};p_{{T}};#varepsilon'.format(key),
                         h2_Numer.GetNbinsX(), h2_Numer.GetXaxis().GetXmin(), h2_Numer.GetXaxis().GetXmax(),
                         h2_Numer.GetNbinsY(), h2_Numer.GetYaxis().GetXmin(), h2_Numer.GetYaxis().GetXmax()))
        h2_Eff[len(h2_Eff)-1].Divide(h2_Numer, h2_Denom)

    print 'Done!!!'
    outFile.Print()

    outFile.Write()
    outFile.Close()

if __name__ == '__main__':

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)

    inFile  = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20131002_220101.root', 'OPEN')

    #bDict = {'bTagEff':('BTruthNumerPt', 'BTruthDenomPt'), 'bMistagEff':('BMistagNumerPt', 'BMistagDenomPt')}
    #ratio_1D(bDict, r.TFile('histos/fcnh_cut1_2012_20130822_183318.root', 'OPEN'), 'inclusive/ttbar')

    eMisQDict = {'leadElectronMisQ':('LeadElecQMisIDNumer', 'LeadElecQMisIDDenom'), 'trailingElectronMisQ':('TrailingElecQMisIDNumer', 'TrailingElecQMisIDDenom')}
    ratio_2D(eMisQDict, 'inclusive/DATA_ELECTRON', inFile)
