#! /usr/bin/env python
import subprocess, shlex, time, math, sys, os, pickle
from array import array
import ROOT as r

seed = int(str(time.time())[6:10])

paramFile   = open('scripts/fcncParams.pkl', 'rb')
scales      = pickle.load(paramFile)
styles      = pickle.load(paramFile)
combos      = pickle.load(paramFile)
removes     = pickle.load(paramFile)
categories  = pickle.load(paramFile)
systematics = pickle.load(paramFile)


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


def card_producer(yields, categories, backgrounds, cardFile):
     
    nChannels   = len(categories)
    nBG         = len(backgrounds)

    nObs   = []
    nSigBG = []
    for category in categories:
        nObs.append(str(yields[category]._yields['data_obs']))

        nSigBG.append('{0:.2f}'.format(yields[category]._yields['fcnh']))
        for bg in backgrounds:
            nSigBG.append('{0:.2f}'.format(yields[category]._yields[bg]))


    ### Need to set observed to nBG ###

    cardFile.write('imax {0}  number of channels\n'.format(nChannels))
    cardFile.write('jmax *  number of backgrounds (\'*\' = automatic)\n')
    cardFile.write('kmax *  number of nuisance parameters (sources of systematical uncertainties)\n')
    cardFile.write('------------\n')
    cardFile.write('bin             {0}\n'.format('\t'.join(categories)))
    cardFile.write('observation     {0}\n'.format('\t'.join(nObs)))
    cardFile.write('------------\n')
    cardFile.write('bin             {0}\n'.format(''.join([(nBG+1)*'{0}\t'.format(category) for category in categories])))
    cardFile.write('process         {0}\n'.format(''.join(nChannels*['{0}\t'.format(i) for i in range(nBG+1)])))
    cardFile.write('process         {0}\n'.format(''.join(nChannels*['{0}\t'.format(sigBG) for sigBG in ['fcnh'] + backgrounds])))
    cardFile.write('rate            {0}\n'.format('\t'.join(nSigBG)))
    cardFile.write('------------\n')

    cardFile.write('lumi    lnN     1.026    1.026  --      --      1.026    1.026  --      --      1.026   1.026   --      --      # 2.6% lumi uncertainty, affects signal and MC-driven background\n')
    cardFile.write('jes     lnN     1.005   1.005   --      --      1.005   1.005   --      --      1.005   1.005   --      --      # 0.5% JES uncertainty (verify this)\n')
    cardFile.write('MET     lnN     1.04    1.04    --      --      1.04    1.04    --      --      1.04    1.04    --      --      # 1% MET resolution uncertainty (verify this)\n')
    cardFile.write('WZ      lnN     --      1.01    --      --      --      1.01    --      --      --      1.01    --      --      # 1% error on irreducible background from WZ contribution\n')
    cardFile.write('ttbar   lnN     1.01    --      --      --      1.01    --      --      --      1.01    --      --      --      # 1% uncertainty on ttbar simulation\n')
    cardFile.write('mu_eff  lnN     1.014   1.014   --      --      1.014   1.014   --      --      1.02    1.02    --      --      # 5% uncertainty on efficiencies (correlated just for simplicity)\n')
    cardFile.write('el_eff  lnN     1.008   1.008   --      --      1.008   1.008   --      --      --      --      --      --      # 5% uncertainty on efficiencies (correlated just for simplicity)\n')
    cardFile.write('qFlips  lnN     --      --      --      1.20     --      --      --      1.20     --      --      --      --    # 1% uncertainty charge flips\n')
    cardFile.write('fakes   lnN     --      --      1.20     --      --      --      1.20     --      --      --      1.20    --    # 5% flat uncertainty on fake rate estimation\n')
    cardFile.write('pileup  lnN     1.01    1.01    --      --      1.01    1.01    --      --      1.01    1.01    --      --      # 1% pileup uncertainty\n')


if __name__ == '__main__':
    if len(sys.argv) > 0:
        batch   = sys.argv[1]
    else:
        print 'A batch and ratio type must be specified.  Otherwise, do some hacking so this thing knows about your inputs.'
        exit()

    r.gStyle.SetOptStat(0)

    doSimpleCuts    = False
    doBinnedScan    = False
    doMask          = True

    lumi        = 19.7e3
    categories  = ['ss_ee', 'ss_emu', 'ss_mumu', 'ss_inclusive']
    backgrounds = ['irr', 'Fakes', 'QFlips']
    datasets    = ['data_obs', 'fcnh'] + backgrounds

    dataDict = {}
    dataDict['data_obs']     = ['DATA_ELECTRON', 'DATA_MUON', 'DATA_MUEG'] # Observed
    dataDict['fcnh']    = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_ZZ_t', 'FCNC_ZZ_tbar', 'FCNC_TauTau_t', 'FCNC_TauTau_tbar'] # signal

    dataDict['irr']     = ['WZJets3LNu', 'ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau', 'ttZ', 'ttW', 'ttG'] # Irreducible backgrounds
    dataDict['Fakes']   = ['muFakes', 'eFakes', 'llFakes'] # Fakes
    dataDict['QFlips']  = ['QFlips'] # electron charge misID

    variable    = 'MetVsHT'
    xBounds     = [60., 300.] # HT
    yBounds     = [10., 140.] # MET

    # input file
    scaleFile   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_{0}.root'.format(batch), 'OPEN')
    histFile    = r.TFile('fcncAnalysis/combined_histos/fcnh_cut3_2012_{0}.root'.format(batch), 'OPEN')
    # Scale factors
    paramFile = open('scripts/fcncParams.pkl', 'rb')

    # prepare output directory
    filePath = 'data/dataCards'
    if not os.path.exists(filePath):
        os.system('mkdir -p '+filePath)
    elif len(os.listdir(filePath)) is not 0:
        os.system('rm -r {0}/*'.format(filePath))


    yields  = {}
    cuts    = []
    for category in categories:
        # Get histograms for different samples 
        histDict = {}
        sumBG = dict(zip(datasets, len(datasets)*[]))
        for dataset in datasets:
            histDict[dataset] = []

            sumBG[dataset] = None
            for i,sample in enumerate(dataDict[dataset]):

                hist = histFile.GetDirectory(category + '/' + sample).Get('h2_' + variable)
                if not hist: continue

                ### Do scaling of MC samples ###
                if dataset not in ['data_obs', 'QFlips', 'Fakes']:
                    yieldHist   = scaleFile.GetDirectory('inclusive/' + sample).Get('h1_YieldByCut')
                    nInit       = yieldHist.GetBinContent(1)
                    hist.Scale(scales['2012'][sample]*lumi/nInit)
                #elif sample == 'eFakes' and category == 'ss_ee':
                #    hist.Scale(2.)

                if sumBG[dataset] == None:
                    sumBG[dataset] = hist.Clone()
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

        if category == 'ss_inclusive': ### Draw optimal cut boundaries on BG and Signal plots
            maskFile = open('data/{0}_mask.txt'.format(variable), 'r')

            canvas = r.TCanvas('canvas', 'canvas', 1000, 500)
            pad1 = r.TPad('pad1', '', 0., 0., 0.5, 1., 0)
            pad1.Draw()
            pad2 = r.TPad('pad2', '', 0.5, 0., 1., 1., 0)
            pad2.Draw()

            cuts = [r.TLine(140., 40., 500., 40.), r.TLine(60., 130., 60., 150.)]
            cuts[0].SetLineColor(r.kRed)
            cuts[0].SetLineWidth(3)
            cuts[1].SetLineColor(r.kRed)
            cuts[1].SetLineWidth(3)

            prevXLow = 140.
            for line in maskFile:
                if line[0] == '#':
                    continue

                line = line.strip('\n')
                #print line

                lineData    = line.split(' ')
                yBin        = int(lineData[0])
                xRangeLow   = int(lineData[1])
                xRangeHigh  = int(lineData[2])
                expLimit    = float(lineData[3])

                newCutX = r.TLine(xRangeLow, yBin, prevXLow, yBin)
                newCutX.SetLineColor(r.kRed)
                newCutX.SetLineWidth(3)

                newCutY = r.TLine(xRangeLow, yBin, xRangeLow, yBin+yBinning)
                newCutY.SetLineColor(r.kRed)
                newCutY.SetLineWidth(3)

                prevXLow = xRangeLow

                cuts.append(newCutX)
                cuts.append(newCutY)

            bgHist = sumBG['irr'].Clone() 
            bgHist.Add(sumBG['Fakes'])
            bgHist.Add(sumBG['QFlips'])
            bgHist.SetTitle('Background')
            bgHist.GetXaxis().SetRangeUser(0.,500.)
            bgHist.GetYaxis().SetRangeUser(0.,150.)

            sumBG['fcnh'].SetTitle('Signal')
            sumBG['fcnh'].GetXaxis().SetRangeUser(0.,500.)
            sumBG['fcnh'].GetYaxis().SetRangeUser(0.,150.)

            pad1.cd()
            bgHist.Draw('colz')
            for cut in cuts:
                cut.Draw('same')

            pad2.cd()
            sumBG['fcnh'].Draw('colz')
            for cut in cuts:
                cut.Draw('same')

            canvas.Print('plots/test.pdf')

        if doSimpleCuts:
            for xCutLow in range(int(math.ceil(xBins/2.))):
                for yCutLow in range(int(math.ceil(yBins/2.))):
                    cut = '{0:d}-{1:d}_{2:d}-{3:d}'.format(int(xCutLow*xBinning), int(xCutHigh*xBinning), int(yCutLow*yBinning), int(yCutHigh*yBinning))
                    cuts.append(cut)

                    if cut in yields.keys():
                        yields[cut][category] = CatData(category, datasets)
                    else:
                        yields[cut] = {category:CatData(category, datasets)}

                    for dataset in datasets:
                        if sumBG[dataset]:
                            numEvents = sumBG[dataset].Integral(xCutLow, xBins, yCutLow, yBins)
                        else:
                            numEvents = 0.

                        yields[cut][category].add_data(dataset, numEvents)

        if doBinnedScan:
            for yBin in range(int(yBounds[0]/yBinning), int(yBounds[1]/yBinning)):
                for xCutLow in range(int(xBounds[0]/xBinning), int(xBounds[1]/xBinning-3)):
                    for xCutHigh in range(xCutLow + 3, int(xBounds[1]/xBinning)): # require that the cut cover at least 5 bins (50 GeV)
                        cut = '{0:d}_{1:d}-{2:d}'.format(int(yBin*yBinning), int(xCutLow*xBinning), int(xCutHigh*xBinning))
                        cuts.append(cut)

                        if cut in yields.keys():
                            yields[cut][category] = CatData(category, datasets)
                        else:
                            yields[cut] = {category:CatData(category, datasets)}

                        for dataset in datasets:
                            if sumBG[dataset]:
                                numEvents = sumBG[dataset].Integral(xCutLow, xCutHigh, yBin, yBin)
                            else:
                                numEvents = 0.
                            yields[cut][category].add_data(dataset, numEvents)

        if doMask:
            maskFile = open('data/{0}_mask.txt'.format(variable), 'r')
            cuts.append('optimal')

            numEvents = dict(zip(datasets, len(datasets)*[0.]))
            yields[category] = CatData(category, datasets)

            for line in maskFile:
                if line[0] == '#':
                    continue

                line = line.strip('\n')
                #print line

                lineData = line.split(' ')
                yBin        = int(int(lineData[0])/yBinning)
                xRangeLow   = int(int(lineData[1])/xBinning)
                xRangeHigh  = int(int(lineData[2])/xBinning)
                expLimit    = float(lineData[3])

                if yBin <= 4: 
                    continue

                print xRangeLow, xCutHigh, yBin

                for dataset in datasets:
                    if sumBG[dataset]:
                        numEvents[dataset] += sumBG[dataset].Integral(xRangeLow, xCutHigh, yBin, yBin)
                    else:
                        numEvents[dataset] += 0.

            for dataset in datasets:
                print dataset, numEvents[dataset]
                yields[category].add_data(dataset, numEvents[dataset])
            print '\n'


    if not doMask:
        cuts = sorted(list(set(cuts)))
        for cut in cuts:
            dataCard = open('{0}/{1}_{2}.txt'.format(filePath, variable, cut), 'w')
            card_producer(yields[cut], categories, backgrounds, dataCard)
    else:
        dataCard = open('{0}/{1}_{2}.txt'.format(filePath, variable, 'optimal'), 'w')
        card_producer(yields, categories, backgrounds, dataCard)




