#! /usr/bin/env python
import subprocess, shlex, time, math, sys, os, pickle
from array import array
import ROOT as r

paramFile   = open('scripts/fcncParams.pkl', 'rb')
scales      = pickle.load(paramFile)
styles      = pickle.load(paramFile)
combos      = pickle.load(paramFile)
removes     = pickle.load(paramFile)
cats        = pickle.load(paramFile)
systematics = pickle.load(paramFile)

class CatData():
    def __init__(self, category, samples):
        self._category  = ''
        self._samples   = samples
        self._yields    = dict(zip(self._samples, len(self._samples)*[0.]))

    def add_data(self, sample, sampleYield):
        self._yields[sample] = sampleYield

class CardProducer():
    def __init__(self, categories, datasets, dataDict, scaleFile):
        self._categories    = categories
        self._datasets      = datasets
        self._dataDict      = dataDict
        self._scaleFile     = scaleFile
        self._lumi          = 1000 # initialize to 1/fb

    def set_luminosity(self, lumi):
        self._lumi = lumi

    def card_producer(self, yields, backgrounds, cardFile):
         
        nChannels   = len(self._categories)
        nBG         = len(backgrounds)

        nObs   = []
        nSigBG = []

        for category in self._categories:
            nObs.append(str(yields[category]._yields['data_obs']))
            nSigBG.append('{0:.2f}'.format(yields[category]._yields['fcnh']))
            for bg in backgrounds:
                nSigBG.append('{0:.2f}'.format(yields[category]._yields[bg]))
            #    nSigBG.append('{0:.2f}'.format(yields[category]._yields[bg]))


        ### Need to set observed to nBG ###

        cardFile.write('imax {0}  number of channels\n'.format(nChannels))
        cardFile.write('jmax *  number of backgrounds (\'*\' = automatic)\n')
        cardFile.write('kmax *  number of nuisance parameters (sources of systematical uncertainties)\n')
        cardFile.write('------------\n')
        cardFile.write('bin             {0}\n'.format('\t'.join(self._categories)))
        cardFile.write('observation     {0}\n'.format('\t'.join(nObs)))
        cardFile.write('------------\n')
        cardFile.write('bin             {0}\n'.format(''.join([(nBG+1)*'{0}\t'.format(category) for category in self._categories])))
        cardFile.write('process         {0}\n'.format(''.join(nChannels*['{0}\t'.format(i) for i in range(nBG+1)])))
        cardFile.write('process         {0}\n'.format(''.join(nChannels*['{0}\t'.format(sigBG) for sigBG in ['fcnh'] + backgrounds])))
        cardFile.write('rate            {0}\n'.format('\t'.join(nSigBG)))
        cardFile.write('------------\n')

        systList = ['lumi', 'jes', 'MET', 'pileup', 'mu_eff', 'el_eff', 'mu_trig', 'el_trig', 'qFlips', 'mufakes', 'elFakes', 'PDF', 'ttbar', 'ttW', 'ttZ', 'WZ']
        systList = ['lumi', 'jes', 'MET', 'pileup', 'mu_eff', 'el_eff', 'mu_trig', 'el_trig', 'mufakes', 'elFakes', 'PDF', 'ttbar', 'ttW', 'ttZ', 'WZ']
        for i, systematic in enumerate(systList):
            cardFile.write('{0}\tlnN\t'.format(systematic))
            for category in self._categories:
                for sample in ['fcnh'] + backgrounds:
                    if sample not in systematics[category].keys(): 
                        continue
                    else:
                        syst = systematics[category][sample][i]
                        if syst != 0.:
                            cardFile.write('{0}\t'.format(1. + systematics[category][sample][i]))
                        else:
                            cardFile.write('--\t')
            cardFile.write('\n')

            #print '{0}\t lnN {1}'.format(systematic, ''.join(['{0}'.format(systematics[category]) for category in categories]))
            #cardFile.write('{0}\t lnN {1}\n'.format('')

    def get_hists_for_yields(self, histFile, variable, category):
        # N.B. Output will be normalized to 1 /fb

        # Get histograms for different samples 
        sumBG = dict(zip(self._datasets, len(self._datasets)*[]))
        for dataset in self._datasets:
            sumBG[dataset] = None
            if dataset in self._dataDict.keys():
                for sample in self._dataDict[dataset]:
                    hist = histFile.GetDirectory(category + '/' + sample).Get(variable)
                    if not hist: continue

                    ### Do scaling of MC samples ###
                    if dataset not in ['data_obs', 'QFlips', 'Fakes']:
                        yieldHist   = self._scaleFile.GetDirectory('inclusive/' + sample).Get('h1_YieldByCut')
                        nInit       = yieldHist.GetBinContent(1)
                        hist.Scale(self._lumi*scales['2012'][sample]/nInit)

                    if sumBG[dataset] == None:
                        sumBG[dataset] = hist.Clone()
                    else:
                        sumBG[dataset].Add(hist)
            else:
                hist = histFile.GetDirectory(category + '/' + dataset).Get(variable)

                if not hist: continue

                if dataset not in ['data_obs', 'QFlips', 'eFakes', 'muFakes', 'llFakes']:
                    yieldHist   = self._scaleFile.GetDirectory('inclusive/' + dataset).Get('h1_YieldByCut')
                    nInit       = yieldHist.GetBinContent(1)
                    hist.Scale(self._lumi*scales['2012'][dataset]/nInit)

                if hist:
                    sumBG[dataset] = hist.Clone()

        return sumBG


if __name__ == '__main__':

    if len(sys.argv) > 0:
        batch   = sys.argv[1]
    else:
        print 'A batch and ratio type must be specified.  Otherwise, do some hacking so this thing knows about your inputs.'
        exit()

    r.gStyle.SetOptStat(0)

    lumi        = 19.7e3
    cutLevel    = 8
    #categories  = ['ss_ee', 'ss_emu', 'ss_mumu']
    categories  = ['3l_eee', '3l_eemu', '3l_emumu', '3l_mumumu']
    backgrounds = ['Irreducible', 'ttW', 'ttZ', 'WZJets3LNu', 'muFakes', 'eFakes', 'llFakes']#, 'QFlips']
    datasets    = ['data_obs', 'fcnh'] + backgrounds

    dataDict = {}
    dataDict['data_obs']    = ['DATA_ELECTRON', 'DATA_MUON', 'DATA_MUEG'] # Observed
    dataDict['fcnh']        = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_ZZ_t', 'FCNC_ZZ_tbar', 'FCNC_TauTau_t', 'FCNC_TauTau_tbar'] # signal
    dataDict['Irreducible'] = ['ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau', 'ttG'] # Irreducible backgrounds
    dataDict['Fakes']       = ['muFakes', 'eFakes', 'llFakes'] # Fakes

    # input file
    scaleFile   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_{0}.root'.format(batch), 'OPEN')
    histFile    = r.TFile('fcncAnalysis/combined_histos/fcnh_cut3_2012_{0}.root'.format(batch), 'OPEN')

    # Scale factors
    paramFile = open('scripts/fcncParams.pkl', 'rb')

    # prepare output directory
    filePath = 'data/dataCards'    
    if not os.path.exists(filePath):
        os.system('mkdir -p '+filePath)

    cardMaker = CardProducer(categories, datasets, dataDict, scaleFile)
    cardMaker.set_luminosity(lumi)

    # get yields for each dataset for the given category
    yields = {}
    for category in categories:
        # Get histograms for different samples 
        sumBG = cardMaker.get_hists_for_yields(histFile, 'h1_YieldByCut', category)
        yields[category] = CatData(category, datasets)
        for dataset in datasets:
            if sumBG[dataset]:
                yields[category].add_data(dataset, sumBG[dataset].GetBinContent(cutLevel))

        print category

    
    dataCard = open('{0}/{1}_{2}.txt'.format(filePath, batch, '3l'), 'w')
    cardMaker.card_producer(yields, backgrounds, dataCard)
