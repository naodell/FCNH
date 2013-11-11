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
    g_Eff.SetTitle('{0};{1};#varepsilon'.format(outName, h1_Numer.GetXaxis().GetTitle(), h1_Numer.GetYaxis().GetTitle()))

    return g_Eff

def make_graph_ratio_2D(outName, h2_Numer, h2_Denom):
    ### Mostly useful for 2D efficiencies, returns an array of
    ### TGraphAsymmErrors

    gList = []

    for i in range(h2_Numer.GetYaxis().GetNbins()):

        gList.append(r.TGraphAsymmErrors())
        gList[i].Divide(h2_Numer.ProjectionX('numer_{0}'.format(i+1), i, i+1), h2_Denom.ProjectionX('denom_{0}'.format(i+1), i, i+1))
        gList[i].SetName('g_{0}_{1}'.format(outName, i+1))
        gList[i].SetTitle('{0}_{1};{2};#varepsilon'.format(outName, i+1, h2_Numer.GetXaxis().GetTitle()))

    return gList


class RatioMaker(AnalysisTools):
    ### For producing efficiency weight files ###
    def __init__(self, inputFileName, outFileName, scale = 1.):
        AnalysisTools.__init__(self, inputFileName, scale)
        self._outFile       = r.TFile(outFileName, 'RECREATE')
        self._ratioDict1D   = {}
        self._ratioDict2D   = {}
        self._hists         = []

    def write_outfile(self):
        self._outFile.Write()
        self._outFile.Close()

    def set_ratio_1D(self, ratioDict):
        self._ratioDict1D = ratioDict

    def set_ratio_2D(self, ratioDict):
        self._ratioDict2D = ratioDict

    def make_category_directory(self):

        if not self._outFile.GetDirectory(self._category):
            self._outFile.mkdir(self._category)
            self._outFile.cd(self._category)

                
    def make_1D_ratios(self, ratioSample, bgSample = '', removePass = True): 
        ### make ratios for all variables specified in ratioDict1D.  Sample
        ### combinations should be specified in combineDict in parameters.py.
        ### bgSample is subtracted off of the inputs for the ratio.

        self.make_category_directory()

        for key,value in self._ratioDict1D.iteritems():
            h1_Numer    = self.combine_samples(value[0], ratioSample) 
            h1_Denom    = self.combine_samples(value[1], ratioSample) 

            if removePass:
                h1_Denom.Add(h1_Numer, -1)

            if bgSample is not '':
                h1_bgNumer  = self.combine_samples(value[0], bgSample) 
                h1_bgDenom  = self.combine_samples(value[1], bgSample) 

                if removePass:
                    h1_bgDenom.Add(h1_bgNumer, -1)

                h1_Numer.Add(h1_bgNumer, -1)
                h1_Denom.Add(h1_bgDenom, -1)

            g_Ratio = make_graph_ratio_1D(key, h1_Numer, h1_Denom)

            self._outFile.GetDirectory(self._category).Add(g_Ratio)


    def make_2D_ratios(self, ratioSample, bgSample = '', doProjections = True, removePass = True): 
        ### make ratios for all variables specified in ratioDict2D.  Sample
        ### combinations should be specified in combineDict in parameters.py.
        ### bgSample is subtracted off of the inputs for the ratio. Produces a
        ### graph for each row in the input 2D histogram

        self.make_category_directory()

        for key,value in self._ratioDict2D.iteritems():
            h2_Numer    = self.combine_samples(value[0], ratioSample, histType = '2D') 
            h2_Denom    = self.combine_samples(value[1], ratioSample, histType = '2D') 

            if removePass:
                h2_Denom.Add(h2_Numer, -1)

            if bgSample is not '':
                h2_bgNumer  = self.combine_samples(value[0], bgSample, histType = '2D') 
                h2_bgDenom  = self.combine_samples(value[1], bgSample, histType = '2D') 

                if removePass:
                    h2_bgDenom.Add(h2_bgNumer, -1)

                h2_Numer.Add(h2_bgNumer, -1.)
                h2_Denom.Add(h2_bgDenom, -1.)

            if doProjections:
                g_RatioList = make_graph_ratio_2D(key, h2_Numer, h2_Denom)

                for g_Ratio in g_RatioList:
                    self._outFile.GetDirectory(self._category).Add(g_Ratio)

            else:
                h2_Eff = r.TH2D('h2_{0}'.format(key), '{0};;'.format(key),
                                 h2_Numer.GetNbinsX(), h2_Numer.GetXaxis().GetXmin(), h2_Numer.GetXaxis().GetXmax(),
                                 h2_Numer.GetNbinsY(), h2_Numer.GetYaxis().GetXmin(), h2_Numer.GetYaxis().GetXmax())

                h2_Eff.Divide(h2_Numer, h2_Denom, 1., 1., 'B')

                self._hists.append(h2_Eff)




if __name__ == '__main__':

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)

    if len(sys.argv) > 1:
        batch = sys.argv[1]
    else:
        print 'A batch must specified.  Otherwise, do some hacking so this thing knows about your inputs.'

    doQFlips    = False
    doFakes     = True

    ### For electron charge misID efficiencies ###
    if doQFlips:
        inFile  = 'fcncAnalysis/combined_histos/fcnh_cut1_2012_{0}.root'.format(batch)
        outFile = 'data/eleQMisID_TEST.root'

        ratioMaker = RatioMaker(inFile, outFile, scale = 19.7)
        ratioMaker.set_category('inclusive')

        eMisQDict = {
            #'LeadElectronMisQ':('LeadElecQMisIDNumer', 'LeadElecQMisIDDenom'),
            #'TrailingElectronMisQ':('TrailingElecQMisIDNumer', 'TrailingElecQMisIDDenom'),
            'DielectronMisQ':('DileptonQMisIDNumer', 'DileptonQMisIDDenom')
            }

        ratioMaker.set_ratio_2D(eMisQDict)
        ratioMaker.make_2D_ratios('DATA', doProjections = False)

        ratioMaker.write_outfile()

    ### For lepton fake rate measurement ###

    if doFakes:
        inFile  = 'fakeEstimator/histos/{0}.root'.format(batch)
        outFile = 'data/fakeRates_TEST.root'

        ratioMaker = RatioMaker(inFile, outFile, scale = 19.7)
        ratioMaker.get_scale_factors(['FAKE_2l', 'FAKE_3l'], corrected = False)

        fakeDict1D = {
            #'ElectronFakeMet':('EleNumerMet', 'EleDenomMet'),
            'ElectronFakePt':('EleNumerPt', 'EleDenomPt'),
            'ElectronFakeEta':('EleNumerEta', 'EleDenomEta'),
            #'MuonFakePt_Even':('MuPassLepPt', 'MuProbeLepPt'),
            #'MuonFakeEta_Even':('MuPassLepEta', 'MuProbeLepEta'),
            #'MuonFakeMet':('MuNumerMet', 'MuDenomMet'),
            'MuonFakePt':('MuNumerPt', 'MuDenomPt'),
            'MuonFakeEta':('MuNumerEta', 'MuDenomEta')
        }

        fakeDict2D = {
            #'MuonFakePt_Even':('MuPassLepPt', 'MuProbeLepPt'),
            #'MuonFakeEta_Even':('MuPassLepEta', 'MuProbeLepEta'),
            'MuonFake':('MuNumer', 'MuDenom'),
            'ElectronFake':('EleNumer', 'EleDenom'),
        }

        fakeCategories = [
                            'QCD2l_inclusive', 'ZPlusJet_inclusive',
                            #'QCD2l_low_met', 'QCD2l_high_met',
                            #'ZPlusJet_low_met', 'ZPlusJet_high_met'
                         ]

        for category in fakeCategories:
            ratioMaker.set_category(category)

            bgType =''
            if category.split('_', 1)[0] == 'QCD2l':
                bgType = 'FAKE_2l'
            elif category.split('_', 1)[0] == 'ZPlusJet':
                bgType = 'FAKE_3l'

            ratioMaker.set_ratio_1D(fakeDict1D)
            ratioMaker.make_1D_ratios('DATA', bgType)

            ratioMaker.set_ratio_2D(fakeDict2D)
            ratioMaker.make_2D_ratios('DATA', bgType)

        ratioMaker.write_outfile()
