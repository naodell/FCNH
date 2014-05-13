#! /usr/bin/env python
from array import array
from AnalysisTools import *
import subprocess
import ROOT as r

def make_graph_ratio_1D(outName, h1_Numer, h1_Denom):
    ### Mostly useful for 1D efficiencies ###

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
        h1_Numer = h2_Numer.ProjectionX('h_Numer_{0}_{1}'.format(outName, i+1), i, i+1)
        h1_Denom = h2_Denom.ProjectionX('h_Denom_{0}_{1}'.format(outName, i+1), i, i+1)

        for j in range(h1_Numer.GetNbinsX()):
            numerContent = h1_Numer.GetBinContent(j+1)
            if numerContent < 0:
                h1_Numer.SetBinContent(i+1, 0.)

            denomContent = h1_Denom.GetBinContent(j+1)
            if denomContent < 0:
                h1_Denom.SetBinContent(i+1, 0.)

        gList.append(r.TGraphAsymmErrors())
        gList[i].Divide(h1_Numer, h1_Denom)
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

                
    def make_1D_ratios(self, ratioSample, bgSample = '', categories = []):
        ### make ratios for all variables specified in ratioDict1D.  Sample
        ### combinations should be specified in combineDict in parameters.py.
        ### bgSample is subtracted off of the inputs for the ratio.

        self.make_category_directory()

        if categories == []:
            categories = [self.get_category(), self.get_category()]

        for key,value in self._ratioDict1D.iteritems():
            self.set_category(categories[0])
            h1_Numer    = self.combine_samples(value[0], ratioSample) 
            self.set_category(categories[1])
            h1_Denom    = self.combine_samples(value[1], ratioSample) 

            if bgSample != '':
                self.set_category(categories[0])
                h1_bgNumer  = self.combine_samples(value[0], bgSample) 
                self.set_category(categories[1])
                h1_bgDenom  = self.combine_samples(value[1], bgSample) 

                h1_Numer.Add(h1_bgNumer, -1.)
                h1_Denom.Add(h1_bgDenom, -1.)

            # Make sure there are no negative weighted histogram bins
            for i in range(h1_Numer.GetNbinsX()):
                binContent = h1_Numer.GetBinContent(i+1),
                if binContent[0] < 0.:
                    h1_Numer.SetBinContent(i+1,0.)

            for i in range(h1_Denom.GetNbinsX()):
                binContent = h1_Denom.GetBinContent(i+1),
                if binContent[0] < 0.:
                    h1_Denom.SetBinContent(i+1,0.)

            ### Save ratios to 2D histograms
            h1_Eff = r.TH1D('h1_{0}'.format(key), '{0};;'.format(key), h1_Numer.GetNbinsX(), h1_Numer.GetXaxis().GetXmin(), h1_Numer.GetXaxis().GetXmax())
            h1_Eff.Divide(h1_Numer, h1_Denom, 1., 1., 'B')
            #h2_Eff.Print("range")
            self._hists.append(h1_Eff)


            g_Ratio = make_graph_ratio_1D(key, h1_Numer, h1_Denom)
            self._outFile.GetDirectory(self._category).Add(g_Ratio)



    def make_2D_ratios(self, ratioSample, bgSample = '', doProjections = True): 
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

            #h2_Numer.Print("range")
            #print ''
            #h2_Denom.Print("range")

            ### Set negative entries to 0
            for binX in range(h2_Numer.GetNbinsX()):
                for binY in range(h2_Numer.GetNbinsY()):
                    if h2_Numer.GetBinContent(binX+1, binY+1) < 0.:
                        h2_Numer.SetBinContent(binX+1, binY+1, 0.)

                    if h2_Denom.GetBinContent(binX+1, binY+1) < 0.:
                        h2_Denom.SetBinContent(binX+1, binY+1, 0.)

            ### Save ratios to 2D histograms
            h2_Eff = r.TH2D('h2_{0}'.format(key), '{0};;'.format(key),
                             h2_Numer.GetNbinsX(), h2_Numer.GetXaxis().GetXmin(), h2_Numer.GetXaxis().GetXmax(),
                             h2_Numer.GetNbinsY(), h2_Numer.GetYaxis().GetXmin(), h2_Numer.GetYaxis().GetXmax())
            h2_Eff.Divide(h2_Numer, h2_Denom, 1., 1., 'B')
            #h2_Eff.Print("range")
            self._hists.append(h2_Eff)

            ### Save ratios to 1D graphs (TGraphAsymmErrors)
            g_RatioList = make_graph_ratio_2D(key, h2_Numer, h2_Denom)
            for g_Ratio in g_RatioList:
                #g_Ratio.Print()
                self._outFile.GetDirectory(self._category).Add(g_Ratio)


    def combine_rates(self, categories):
        ### Combines fake rates for differenct categories

        for category in categories:
            sourceDir = self._outFile.GetDirectory(categories[0])


    def charge_flip_fitter(self, ratioSample, nToys = 10):
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

            prob0 = [[0.5*h2_Eff.GetBinContent(i+1, i+1) for i in range(nBinsX)], [0.5*h2_Eff.GetBinError(i+1, i+1) for i in range(nBinsX)]]

            for toy in range(nToys):

                ### Now get possible values of p(i) from p(i) = P(i,j) - p(j) and reiterate
                probs = [[0. for i in range(nBinsX)], [0. for i in range(nBinsX)]]

                for binX in range(nBinsX):
                    for binY in range(nBinsY):
                        #if binX == binY: continue

                        binContentXY    = h2_Eff.GetBinContent(binX+1, binY+1) 
                        binErrorXY      = h2_Eff.GetBinError(binX+1, binY+1) 
                        errX            = binErrorXY*binErrorXY + prob0[1][binY]*prob0[1][binY]
                        errY            = binErrorXY*binErrorXY + prob0[1][binX]*prob0[1][binX]
                        if binContentXY != 0:
                            if prob0[0][binY] != 0:
                                #print binContentXY, prob0[binY][0], errX

                                probs[0][binX] += (binContentXY - prob0[0][binY])/errX
                                probs[1][binX] += 1./errX
                            if prob0[0][binX] != 0:
                                probs[0][binY] += (binContentXY - prob0[0][binX])/errY
                                probs[1][binY] += 1./errY

                for binX in range(nBinsX):
                    #print probs[binX], probs[0][binX]/nBinsX
                    if probs[1][binX] != 0:
                        prob0[0][binX] = probs[0][binX]/probs[1][binX]
                        prob0[1][binX] = 1./sqrt(probs[1][binX])


            #ptBins = [15., 25., 37.5, 57.5, 85., 125.]
            ptBins = [15., 27.5, 55., 115.]
            g_ProbB = r.TGraphErrors(len(ptBins), array('f', ptBins),  array('f', prob0[0][:nBinsX/2]), array('f', [0.1 for bin in ptBins]), array('f', prob0[1][:nBinsX/2]))
            g_ProbB.SetName('g_QFlipB')
            g_ProbB.SetTitle('barrel electron charge flips;iPt;#varepsilon')

            g_ProbE = r.TGraphErrors(len(ptBins), array('f', ptBins),  array('f', prob0[0][nBinsX/2:]), array('f', [0.1 for bin in ptBins]), array('f', prob0[1][nBinsX/2:]))
            g_ProbE.SetName('g_QFlipE')
            g_ProbE.SetTitle('endcap electron charge flips;iPt;#varepsilon')

            self._outFile.GetDirectory(self._category).Add(g_ProbB)
            self._outFile.GetDirectory(self._category).Add(g_ProbE)


if __name__ == '__main__':

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)

    if len(sys.argv) > 1:
        type    = sys.argv[1]
        batch   = sys.argv[2]
    else:
        print 'A batch and ratio type must be specified.  Otherwise, do some hacking so this thing knows about your inputs.'
        exit()

    ### For AIC rates ###
    if type == 'AIC':
        inFile  = 'fcncAnalysis/combined_histos/fcnh_cut1_2012_{0}.root'.format(batch)
        outFile = 'data/AIC_TEST.root'

        ratioMaker = RatioMaker(inFile, outFile, scale = 19.7)
        ratioMaker.set_category('inclusive')
        ratioMaker.get_scale_factors(['AIC_BG'], corrected = False)

        ratioMaker.set_ratio_1D({'mumumu':('TrileptonMass_AIC', 'DimuonPhotonMass_AIC')})
        ratioMaker.make_1D_ratios('DATA', bgSample = '', categories = ['3l_mumumu', 'inclusive'])
        ratioMaker.set_ratio_1D({'emumu':('TrileptonMass_AIC', 'DimuonPhotonMass_AIC')})
        ratioMaker.make_1D_ratios('DATA', bgSample = '', categories = ['3l_emumu', 'inclusive'])

        ratioMaker.set_ratio_1D({'eemu':('TrileptonMass_AIC', 'DimuonPhotonMass_AIC')})
        ratioMaker.make_1D_ratios('DATA', bgSample = '', categories = ['3l_eemu', 'inclusive'])
        ratioMaker.set_ratio_1D({'eee':('TrileptonMass_AIC', 'DimuonPhotonMass_AIC')})
        ratioMaker.make_1D_ratios('DATA', bgSample = '', categories = ['3l_eee', 'inclusive'])

        ratioMaker.write_outfile()

    ### For electron charge misID efficiencies ###
    if type == 'QFlips':
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
        ratioMaker.charge_flip_fitter('DATA_ELECTRON', nToys = 50)
        #ratioMaker.make_2D_ratios('DATA', doProjections = False)

        ratioMaker.write_outfile()

    ### For lepton fake rate measurement ###

    if type == 'Fakes':
        inFile  = 'fakeEstimator/histos/{0}.root'.format(batch)
        outFile = 'data/fakeRates_TEST.root'

        ratioMaker = RatioMaker(inFile, outFile, scale = 19.7)
        #ratioMaker.get_scale_factors(['PROMPT'], corrected = False)

        fakeDict1D = {
            'MuonFakeMet':('MuNumerMet', 'MuDenomMet'),
            'MuonFakePt':('MuNumerPt', 'MuDenomPt'),
            'MuonFakeEta':('MuNumerEta', 'MuDenomEta'),
            #'ElectronFakeMet':('EleNumerMet', 'EleDenomMet'),
            #'ElectronFakePt':('EleNumerPt', 'EleDenomPt'),
            #'ElectronFakeEta':('EleNumerEta', 'EleDenomEta')
        }

        fakeDict2D = {
            'MuonFake':('MuNumer', 'MuDenom'),
            #'ElectronFake':('EleNumer', 'EleDenom')
        }

        fakeCategories = ['QCD2l']#, 'ZPlusJet', 'AntiIso3l']

        for category in fakeCategories:
            ratioMaker.set_category(category)

            #bgType ='PROMPT'
            bgType =''

            ratioMaker.set_ratio_1D(fakeDict1D)
            ratioMaker.make_1D_ratios('DATA', bgType)

            ratioMaker.set_ratio_2D(fakeDict2D)
            ratioMaker.make_2D_ratios('DATA', bgType, doProjections = True)

        ratioMaker.write_outfile()

        # Combined fakeCategory rates
        #ratioMaker.combine_rates(fakeCategories)

