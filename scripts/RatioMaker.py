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
        gList[i].Divide(h2_Numer.ProjectionX('h_Numer_{0}_{1}'.format(outName, i+1), i, i+1), h2_Denom.ProjectionX('h_Denom_{0}_{1}'.format(outName, i+1), i, i+1))
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

    def write_outfile(self, hists = []):
        
        for hist in hists:
            self._outFile.Add(hist)

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


    def make_2D_ratios(self, ratioSample, bgSample = '', doProjections = True, removePass = False): 
        ### make ratios for all variables specified in ratioDict2D.  Sample
        ### combinations should be specified in combineDict in parameters.py.
        ### bgSample is subtracted off of the inputs for the ratio. Produces a
        ### graph for each row in the input 2D histogram

        self.make_category_directory()

        for key,value in self._ratioDict2D.iteritems():
            h2_Numer    = self.combine_samples(value[0], ratioSample, histType = '2D') 
            h2_Denom    = self.combine_samples(value[1], ratioSample, histType = '2D') 

            #print h2_Numer.Integral(),h2_Denom.Integral(), 

            if bgSample != '':
                h2_bgNumer  = self.combine_samples(value[0], bgSample, histType = '2D') 
                h2_bgDenom  = self.combine_samples(value[1], bgSample, histType = '2D') 

                h2_Numer.Add(h2_bgNumer, -1.)
                h2_Denom.Add(h2_bgDenom, -1.)

            ### Set negative entries to 0
            for binX in range(h2_Numer.GetNbinsX()):
                for binY in range(h2_Numer.GetNbinsY()):
                    if h2_Numer.GetBinContent(binX+1, binY+1) < 0.:
                        h2_Numer.SetBinContent(binX+1, binY+1, 0.)

                    if h2_Denom.GetBinContent(binX+1, binY+1) < 0.:
                        h2_Denom.SetBinContent(binX+1, binY+1, 0.)

            if removePass:
                h2_Denom.Add(h2_Numer, -1)

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

    def charge_flip_fitter(self, ratioSample, nToys = 5):
        ### Converts 2D (double electron) charge flip probabilities to 1D
        ### (single electron) charge flip probabilities.  Assumes mapping is
        ### P(i,j) = p(i) + p(j).

        self.make_category_directory()

        for key,value in self._ratioDict2D.iteritems():
            h2_Numer    = self.combine_samples(value[0], ratioSample, histType = '2D') 
            h2_Denom    = self.combine_samples(value[1], ratioSample, histType = '2D') 

            h2_Eff = r.TH2D('h2_{0}'.format(key), '{0};;'.format(key),
                             h2_Numer.GetNbinsX(), h2_Numer.GetXaxis().GetXmin(), h2_Numer.GetXaxis().GetXmax(),
                             h2_Numer.GetNbinsY(), h2_Numer.GetYaxis().GetXmin(), h2_Numer.GetYaxis().GetXmax())

            h2_Eff.Divide(h2_Numer, h2_Denom, 1., 1., 'B')
            self._hists.append(h2_Eff)

            ### Guess at values of p(i) from diagonal bins, i.e., p(i) = 0.5*P(i,j)

            nBinsX = h2_Eff.GetNbinsX()
            nBinsY = h2_Eff.GetNbinsY()

            prob0 = [0.5*h2_Eff.GetBinContent(i+1, i+1) for i in range(nBinsX)]

            for toy in range(nToys):

                ### Now get possible values of p(i) from p(i) = P(i,j) - p(j) and reiterate
                probs = [[] for i in range(nBinsX)]
                for binX in range(nBinsX):
                    for binY in range(nBinsY):
                        if binX == binY: continue

                        binContentXY = h2_Eff.GetBinContent(binX+1, binY+1) 
                        binContentYX = h2_Eff.GetBinContent(binY+1, binX+1) 
                        if binContentXY != 0:
                            if prob0[binY] != 0:
                                probs[binX].append(binContentXY - prob0[binY])
                            if prob0[binX] != 0:
                                probs[binY].append(binContentXY - prob0[binX])
                        if binContentYX != 0:
                            if prob0[binX] != 0:
                                probs[binY].append(binContentYX - prob0[binX])
                            if prob0[binX] != 0: 
                                probs[binX].append(binContentYX - prob0[binY])

                for binX in range(nBinsX):
                    #print probs[binX], sum(probs[binX])/nBinsX
                    prob0[binX] = abs(sum(probs[binX])/len(probs[binX]))

            ptBins = [15., 25., 37.5, 57.5, 85., 125.]
            g_ProbB = r.TGraph(len(ptBins), array('f', ptBins), array('f', prob0[:nBinsX/2]))
            g_ProbB.SetName('g_QFlipB')
            g_ProbB.SetTitle('barrel electron charge flips;iPt;#varepsilon')

            g_ProbE = r.TGraph(len(ptBins), array('f', ptBins), array('f', prob0[nBinsX/2:]))
            g_ProbE.SetName('g_QFlipE')
            g_ProbE.SetTitle('endcap electron charge flips;iPt;#varepsilon')

            self._outFile.GetDirectory(self._category).Add(g_ProbB)
            self._outFile.GetDirectory(self._category).Add(g_ProbE)


if __name__ == '__main__':

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)

    if len(sys.argv) > 1:
        batch = sys.argv[1]
    else:
        print 'A batch must specified.  Otherwise, do some hacking so this thing knows about your inputs.'
        exit()

    doQFlips    = False
    doFakes     = True
    doMetFake   = False

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
        ratioMaker.charge_flip_fitter('DATA')
        #ratioMaker.make_2D_ratios('DATA', doProjections = False)

        ratioMaker.write_outfile()

    ### For lepton fake rate measurement ###

    if doFakes:
        inFile  = 'fakeEstimator/histos/{0}.root'.format(batch)
        outFile = 'data/fakeRates_TEST.root'

        ratioMaker = RatioMaker(inFile, outFile, scale = 19.7)
        ratioMaker.get_scale_factors(['FAKES_2l', 'FAKES_3l'], corrected = False)

        fakeDict1D = {
            #'MuonFakePt_Even':('MuPassLepPt', 'MuProbeLepPt'),
            #'MuonFakeEta_Even':('MuPassLepEta', 'MuProbeLepEta'),
            'MuonFakeMet':('MuNumerMet', 'MuDenomMet'),
            'MuonFakePt':('MuNumerPt', 'MuDenomPt'),
            'MuonFakeEta':('MuNumerEta', 'MuDenomEta'),
            'ElectronFakeMet':('EleNumerMet', 'EleDenomMet'),
            'ElectronFakePt':('EleNumerPt', 'EleDenomPt'),
            'ElectronFakeEta':('EleNumerEta', 'EleDenomEta')
        }

        fakeDict2D = {
            #'MuonFakePt_Even':('MuPassLepPt', 'MuProbeLepPt'),
            #'MuonFakeEta_Even':('MuPassLepEta', 'MuProbeLepEta'),
            'MuonFake':('MuNumer', 'MuDenom'),
            'ElectronFake':('EleNumer', 'EleDenom')
        }

        fakeCategories = ['QCD2l_inclusive', 'ZPlusJet_inclusive']

        for category in fakeCategories:
            ratioMaker.set_category(category)

            bgType =''
            if category.split('_', 1)[0] == 'QCD2l':
                bgType = 'FAKES_2l'
            elif category.split('_', 1)[0] == 'ZPlusJet':
                bgType = 'FAKES_3l'

            ratioMaker.set_ratio_1D(fakeDict1D)
            ratioMaker.make_1D_ratios('DATA', bgType)

            ratioMaker.set_ratio_2D(fakeDict2D)
            ratioMaker.make_2D_ratios('DATA', bgType, doProjections = True, removePass = False)

        # ttH fake rate estimation using low/high met categories

        if doMetFake:
            ratioMaker.set_category('TEST')
            ratioMaker.make_category_directory()

            hists = []
            for key,value in fakeDict2D.iteritems():
                ratioMaker.set_category('QCD2l_low_met')
                h2_Numer_lm    = ratioMaker.combine_samples(value[0], 'DATA', histType = '2D') 
                h2_Denom_lm    = ratioMaker.combine_samples(value[1], 'DATA', histType = '2D') 
                h2_bgNumer_lm  = ratioMaker.combine_samples(value[0], 'FAKES_2l', histType = '2D') 
                h2_bgDenom_lm  = ratioMaker.combine_samples(value[1], 'FAKES_2l', histType = '2D') 

                ratioMaker.set_category('QCD2l_high_met')
                h2_Numer_hm    = ratioMaker.combine_samples(value[0], 'DATA', histType = '2D') 
                h2_Denom_hm    = ratioMaker.combine_samples(value[1], 'DATA', histType = '2D') 
                h2_bgNumer_hm  = ratioMaker.combine_samples(value[0], 'FAKES_2l', histType = '2D') 
                h2_bgDenom_hm  = ratioMaker.combine_samples(value[1], 'FAKES_2l', histType = '2D') 

                h2_bgNumer = h2_bgNumer_lm.Clone()
                h2_bgNumer.Add(h2_bgNumer_hm)

                nBinsX, xMin, xMax, nBinsY, yMin, yMax = \
                                                         h2_Numer_lm.GetNbinsX(), h2_Numer_lm.GetXaxis().GetXmin(), h2_Numer_lm.GetXaxis().GetXmax(), \
                                                         h2_Numer_lm.GetNbinsY(), h2_Numer_lm.GetYaxis().GetXmin(), h2_Numer_lm.GetYaxis().GetXmax()

                h2_fakeRate_lm  = r.TH2D('h2_fakeRate_lm_{0}'.format(key), '{0};;'.format(key), nBinsX, xMin, xMax, nBinsY, yMin, yMax)
                h2_fakeRate_hm  = r.TH2D('h2_fakeRate_hm_{0}'.format(key), '{0};;'.format(key), nBinsX, xMin, xMax, nBinsY, yMin, yMax)
                h2_fakeRate_qcd = r.TH2D('h2_fakeRate_qcd_{0}'.format(key), '{0};;'.format(key), nBinsX, xMin, xMax, nBinsY, yMin, yMax)
                h2_ratioPrompt  = r.TH2D('h2_ratioPrompt_lm_{0}'.format(key), '{0};;'.format(key), nBinsX, xMin, xMax, nBinsY, yMin, yMax)
                h2_Numerator    = r.TH2D('h2_Numerator_lm_{0}'.format(key), '{0};;'.format(key), nBinsX, xMin, xMax, nBinsY, yMin, yMax)
                h2_Denominator  = r.TH2D('h2_Denominator_lm_{0}'.format(key), '{0};;'.format(key), nBinsX, xMin, xMax, nBinsY, yMin, yMax)
                h2_AllOnes      = r.TH2D('h2_AllOnes_lm_{0}'.format(key), '{0};;'.format(key), nBinsX, xMin, xMax, nBinsY, yMin, yMax)

                for n in range(nBinsX):
                    for m in range(nBinsY):
                        h2_AllOnes.SetBinContent(n+1, m+1, 1);
                        h2_AllOnes.SetBinError(n+1, m+1, 0);

                h2_fakeRate_lm.Divide(h2_Numer_lm, h2_Denom_lm, 1., 1., 'B')
                h2_fakeRate_hm.Divide(h2_Numer_hm, h2_Denom_lm, 1., 1., 'B')

                h2_Denom_lm.Multiply(h2_bgNumer_hm)
                h2_Denom_hm.Multiply(h2_bgNumer_lm)
                h2_ratioPrompt.Divide(h2_Denom_lm, h2_Denom_hm, 1., 1.)

                h2_fakeRate_hm.Multiply(h2_ratioPrompt)
                h2_Numerator.Add(h2_fakeRate_lm, h2_fakeRate_hm, 1., -1.)
                h2_Denominator.Add(h2_AllOnes, h2_ratioPrompt, 1., -1.)

                h2_fakeRate_qcd.Divide(h2_Numerator, h2_Denominator, 1., 1., 'B')

                hists.append(h2_fakeRate_qcd)

        ratioMaker.write_outfile()
