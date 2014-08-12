#! /usr/bin/env python
import subprocess, shlex, time, math, sys, os, pickle
from array import array
import ROOT as r

seed = int(str(time.time())[6:10])

'''
Carries out a random grid search to optimize rectangular cuts in 2D plane.
'''

class CatData():
    def __init__(self, category, samples):
        self._category  = ''
        self._samples   = samples
        self._yields    = dict(zip(self._samples, len(self._samples)*[0.]))

    def add_data(self, sample, sampleYield):
        self._yields[sample] = sampleYield


def card_producer(yields, categories, outFile):
     
    nChannels   = 1

    ### Need to set observed to nBG ###

    outFile.write('imax {0}  number of channels'.format(nChannels))
    outFile.write('jmax *  number of backgrounds (\'*\' = automatic)')
    outFile.write('kmax *  number of nuisance parameters (sources of systematical uncertainties)')
    outFile.write('------------')
    outFile.write('bin             {0}'.format(categories))
    outFile.write('observation     {0}'.format(nObs))
    outFile.write('------------')
    outFile.write('bin             {0}'.format(categories))
    outFile.write('process         {0}'.format(backgrounds*nCats))
    outFile.write('process         {0}'.format('0.\t'*nCats))
    outFile.write('rate            {0}'.format(nBG))
    outFile.write('------------')

    #lumi    lnN     1.026    1.026  --      --      1.026    1.026  --      --      1.026   1.026   --      --      # 2.6% lumi uncertainty, affects signal and MC-driven background
    #jes     lnN     1.005   1.005   --      --      1.005   1.005   --      --      1.005   1.005   --      --      # 0.5% JES uncertainty (verify this)
    #MET     lnN     1.04    1.04    --      --      1.04    1.04    --      --      1.04    1.04    --      --      # 1% MET resolution uncertainty (verify this)
    #WZ      lnN     --      1.01    --      --      --      1.01    --      --      --      1.01    --      --      # 1% error on irreducible background from WZ contribution
    #ttbar   lnN     1.01    --      --      --      1.01    --      --      --      1.01    --      --      --      # 1% uncertainty on ttbar simulation
    #mu_eff  lnN     1.014   1.014   --      --      1.014   1.014   --      --      1.014   1.014   --      --      # 5% uncertainty on efficiencies (correlated just for simplicity)
    #el_eff  lnN     1.008   1.008   --      --      1.008   1.008   --      --      1.008   1.008   --      --      # 5% uncertainty on efficiencies (correlated just for simplicity)
    #qFlips  lnN     --      --      --      1.20     --      --      --      1.20     --      --      --      --    # 1% uncertainty charge flips
    #fakes   lnN     --      --      1.20     --      --      --      1.20     --      --      --      1.20    --    # 5% flat uncertainty on fake rate estimation
    #pileup  lnN     1.01    1.01    --      --      1.01    1.01    --      --      1.01    1.01    --      --      # 1% pileup uncertainty


if __name__ == '__main__':
    if len(sys.argv) > 0:
        batch   = sys.argv[1]
    else:
        print 'A batch and ratio type must be specified.  Otherwise, do some hacking so this thing knows about your inputs.'
        exit()

    categories  = ['ss_ee', 'ss_emu', 'ss_mumu']
    datasets    = ['obs', 'fcnh', 'irr']#, 'Fakes', 'QFlips']
    yields  = {}

    dataDict = {}
    dataDict['obs']     = ['DATA_ELECTRON', 'DATA_MUON', 'DATA_MUEG'] # Observed
    dataDict['fcnh']    = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_ZZ_t', 'FCNC_ZZ_tbar', 'FCNC_TauTau_t', 'FCNC_TauTau_tbar'] # signal
    dataDict['irr']     = ['WZJets3LNu', 'ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau', 'ttZ', 'ttW', 'ttG'] # Irreducible backgrounds
    dataDict['Fakes']   = ['muFakes', 'eFakes', 'llFakes'] # Fakes
    dataDict['QFlips']  = ['QFlips'] # electron charge misID

    variable = 'metVsHt'

    # input file
    hFile = r.TFile('fcncAnalysis/combined_histos/fcnh_cut3_2012_{0}.root'.format(batch), 'OPEN')
    # Scale factors
    paramFile = open('scripts/fcncParams.pkl', 'rb')
    scales    = pickle.load(paramFile)

    for category in categories:
        # Get histograms for different samples 
        histDict = {}
        sumBG = dict(zip(datasets, len(datasets)*[]))
        for dataset in datasets:
            histDict[dataset] = []

            sumBG[dataset] = r.TH2D()
            for i,sample in enumerate(dataDict[dataset]):
                #print category, sample, variable
                hist = hFile.GetDirectory(category + '/' + sample).Get('h2_' + variable)

                if not hist:
                    continue

                #print sample
                hist.Scale(scales['2012'][sample])

                if i is 0:
                    sumBG[dataset] = hist
                else:
                    sumBG[dataset].Add(hist)

        xLow    = sumBG[datasets[0]].GetXaxis().GetXmin() 
        xHigh   = sumBG[datasets[0]].GetXaxis().GetXmax() 
        yLow    = sumBG[datasets[0]].GetYaxis().GetXmin() 
        yHigh   = sumBG[datasets[0]].GetYaxis().GetXmax() 
        xBins   = sumBG[datasets[0]].GetNbinsX()
        yBins   = sumBG[datasets[0]].GetNbinsY()

        xBinning = (xHigh - xLow)/xBins
        yBinning = (yHigh - yLow)/yBins
        xCutHigh = xBins
        yCutHigh = yBins

        for xCutLow in range(1, int(math.ceil(xBins/2.))):
            for yCutLow in range(1,int(math.ceil(yBins/2.))):
                cut = '{0:d}-{1:d}_{2:d}-{3:d}'.format(int(xCutLow*xBinning), int(xCutHigh*xBinning), int(yCutLow*yBinning), int(yCutHigh*yBinning))

                if cut in yields.keys():
                    yields[cut][category] = CatData(category, datasets)
                else:
                    yields[cut] = {category:CatData(category, datasets)}

                for dataset in datasets:
                    numEvents = sumBG[dataset].Integral(xCutLow, xBins, yCutLow, yBins)
                    yields[cut][category].add_data(dataset, numEvents)

    print yields['480-1000_170-350']['ss_ee']._yields

    #for i in range(nCuts):
    #    dataCard = open('test.txt', 'w')
    #    card_producer(yields, categories, dataCard)

