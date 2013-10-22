#! /usr/bin/env python
from array import array
from AnalysisTools import *
import subprocess
import ROOT as r

def make_graph_ratio_1D(outName, h1_Numer, h1_Denom):
    ### Mostly useful for 1D efficiencies 

    g_Eff = r.TGraphAsymmErrors()
    g_Eff.Divide(h1_Numer, h1_Denom)
    g_Eff.SetName('g_{0}'.format(outName))
    g_Eff.SetTitle('{0};{1};{2}'.format(outName, h1_Numer.GetXaxis().GetTitle(), h1_Numer.GetYaxis().GetTitle()))

    return g_Eff

def make_graph_ratio_2D(outName, h2_Numer, h2_Denom):
    ### Mostly useful for 2D efficiencies, returns an array of
    ### TGraphAsymmErrors

    gList = []

    for i in range(h2_Numer.GetXaxis().GetNbins()):
        g_Eff.apppend(r.TGraphAsymmErrors())
        g_Eff[i].Divide(h1_Numer, h1_Denom)
        g_Eff[i].SetName('g_{0}'.format(key))
        g_Eff[i].SetTitle('{0};{1};{2}'.format(outName, h1_Numer.GetXaxis().GetTitle(), h1_Numer.GetYaxis().GetTitle()))

    return gList


class RatioMaker(AnalysisTools):
    ### For producing efficiency weight files ###
    def __init__(self, inputFileName, outFileName, scale = 1.):
        AnalysisTools.__init__(self, inputFileName, scale)
        self._outFile       = r.TFile(outFileName, 'RECREATE')
        self._ratioDict1D   = {}
        self._ratioDict2D   = {}

    def set_ratio_1D(self, ratioDict):
        self._ratioDict1D = ratioDict

    def set_ratio_2D(self, ratioDict):
        self._ratioDict2D = ratioDict

    def make_1D_ratios(self, ratioSample, bgSample = ''): 
        ### make ratios for all variables specified in ratioDict1D.  Sample
        ### combinations should be specified in combineDict in parameters.py.
        ### bgSample is subtracted off of the inputs for the ratio.

        canvas = r.TCanvas('canvas', 'canvas', 600, 800)

        for key,value in self._ratioDict1D.iteritems():
            h1_Numer    = self.combine_samples(value[0], ratioSample) 
            h1_Denom    = self.combine_samples(value[1], ratioSample) 

            if bgSample is not '':
                h1_bgNumer  = self.combine_samples(value[0], bgSample) 
                h1_bgDenom  = self.combine_samples(value[1], bgSample) 
                h1_Numer.Add(h1_bgNumer, -1)
                h1_Denom.Add(h1_bgDenom, -1)

                g_Ratio = make_graph_ratio_1D(key, h1_Numer, h1_Denom)

            else:
                g_Ratio = make_graph_ratio_1D(key, h1_Numer, h1_Denom)

            self._outFile.Add(g_Ratio)

    def write_outfile(self):
        self._outFile.Write()
        self._outFile.Close()


if __name__ == '__main__':

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)

    if len(sys.argv) > 1:
        batch = sys.argv[1]
    else:
        batch = '20131021_120020'

    ### For electron charge misID efficiencies ###
    #inFile  = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20131003_171426.root', 'OPEN')
    #outFile = r.TFile('data/electronQMisID.root', 'RECREATE')

    #eMisQDict = {
    #    'LeadElectronMisQ':('LeadElecQMisIDNumer', 'LeadElecQMisIDDenom'),
    #    'TrailingElectronMisQ':('TrailingElecQMisIDNumer', 'TrailingElecQMisIDDenom'),
    #    'DielectronMisQ':('DileptonQMisIDNumer', 'DileptonQMisIDDenom')}

    #ratio_2D(eMisQDict, 'inclusive/DATA_ELECTRON', inFile)

    ### For lepton fake rate measurement ###
    inFile  = 'fakeEstimator/histos/{0}.root'.format(batch)
    outFile = 'data/fakeRates_TEST.root'

    ratioMaker = RatioMaker(inFile, outFile)
    ratioMaker.set_category('inclusive')
    ratioMaker.get_scale_factors(['FAKE_BG'], corrected = False)

    fakeDict1D = {
        #'MuonFakePt_Even':('MuPassLepPt', 'MuProbeLepPt'),
        #'MuonFakeEta_Even':('MuPassLepEta', 'MuProbeLepEta'),
        'MuonFakePt':('MuNumerPt', 'MuDenomPt'),
        'MuonFakeEta':('MuNumerEta', 'MuDenomEta')
    }

    ratioMaker.set_ratio_1D(fakeDict1D)
    ratioMaker.make_1D_ratios('DATA', 'FAKE_BG')
    ratioMaker.write_outfile()
