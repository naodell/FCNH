#! /usr/bin/env python
import subprocess, shlex, time, math, sys, os, pickle
from array import array
import ROOT as r
from CardProducer import *

if __name__ == '__main__':
    if len(sys.argv) > 0:
        batch   = sys.argv[1]
    else:
        print 'A batch and ratio type must be specified.  Otherwise, do some hacking so this thing knows about your inputs.'
        exit()

    r.gStyle.SetOptStat(0)
    r.gROOT.SetBatch()

    doSimpleCuts    = False
    doBinnedScan    = False
    doMask          = True

    lumi        = 19.7e3
    #categories  = ['ss_ee', 'ss_emu', 'ss_mumu']
    categories  = ['ss_inclusive']
    backgrounds = ['Irreducible', 'ttW', 'ttZ', 'WZJets3LNu', 'muFakes', 'eFakes', 'llFakes', 'QFlips']
    datasets    = ['data_obs', 'fcnh'] + backgrounds

    dataDict = {}
    dataDict['data_obs']    = ['DATA_ELECTRON', 'DATA_MUON', 'DATA_MUEG'] # Observed
    dataDict['fcnh']        = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_ZZ_t', 'FCNC_ZZ_tbar', 'FCNC_TauTau_t', 'FCNC_TauTau_tbar'] # signal
    dataDict['Irreducible'] = ['ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau', 'ttG'] # Irreducible backgrounds
    dataDict['Fakes']       = ['muFakes', 'eFakes', 'llFakes'] # Fakes

    variable    = 'MetVsHT'
    xBounds     = [60., 300.] # HT
    yBounds     = [10., 140.] # MET

    # input file
    scaleFile   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_{0}.root'.format(batch), 'OPEN')
    histFile    = r.TFile('fcncAnalysis/combined_histos/fcnh_cut3_2012_{0}.root'.format(batch), 'OPEN')
    # Scale factors
    paramFile = open('scripts/fcncParams.pkl', 'rb')

    # prepare output directory
    filePath = 'data/dataCards'.format(batch)
    if not os.path.exists(filePath):
        os.system('mkdir -p '+filePath)
    #elif len(os.listdir(filePath)) is not 0:
    #    os.system('rm -r {0}/*'.format(filePath))

    cardMaker = CardProducer(categories, datasets, dataDict, scaleFile)
    cardMaker.set_luminosity(lumi)

    yields  = {}
    cuts    = []
    for category in categories:
        # Get histograms for different samples 
        sumBG = cardMaker.get_hists_for_yields(histFile, 'h2_' + variable, category)

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

            bgHist = sumBG['Irreducible'].Clone() 
            bgHist.Add(sumBG['muFakes'])
            bgHist.Add(sumBG['eFakes'])
            bgHist.Add(sumBG['llFakes'])
            bgHist.Add(sumBG['QFlips'])
            bgHist.SetTitle('Background')
            bgHist.GetXaxis().SetTitle('H_{T} [GeV]')
            bgHist.GetYaxis().SetTitle('MET [GeV]')
            bgHist.GetXaxis().SetRangeUser(0.,500.)
            bgHist.GetYaxis().SetRangeUser(0.,150.)
            bgHist.GetYaxis().SetTitleOffset(1.3)

            sumBG['fcnh'].SetTitle('Signal')
            sumBG['fcnh'].GetXaxis().SetTitle('H_{T} [GeV]')
            sumBG['fcnh'].GetYaxis().SetTitle('MET [GeV]')
            sumBG['fcnh'].GetXaxis().SetRangeUser(0.,500.)
            sumBG['fcnh'].GetYaxis().SetRangeUser(0.,150.)
            sumBG['fcnh'].GetYaxis().SetTitleOffset(1.3)

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

            binnedYields = [{'DATA':0., 'BG':0., 'fcnh':0.}, {'DATA':0., 'BG':0., 'fcnh':0.}, {'DATA':0., 'BG':0., 'fcnh':0.}] # dataset, bin, nEvents, error

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

                if yBin <= 4: continue

                #print xRangeLow, xCutHigh, yBin 
                if yBin < 7:
                    bin = 0
                elif yBin < 9:
                    bin = 1
                else:
                    bin = 2

                for dataset in datasets:
                    integral    = 0.
                    error       = r.Double()
                    if sumBG[dataset]:
                        integral = sumBG[dataset].IntegralAndError(xRangeLow, xCutHigh, yBin, yBin, error)

                    numEvents[dataset] += 0

                    #binnedYields[yBin-1][dataset] = (numEvents[dataset], error)

                    if dataset == 'data_obs':
                        binnedYields[bin]['DATA'] += integral
                    elif dataset == 'fcnh':
                        binnedYields[bin]['fcnh'] += integral
                    else:
                        binnedYields[bin]['BG']   += integral

            for dataset in datasets:
                yields[category].add_data(dataset, numEvents[dataset])

            for yields in binnedYields:
                for dataset, nEvents in yields.iteritems():
                    print '& {0:.1f} '.format(nEvents),

                print ' \\\\'

    print yields


    if not doMask:
        cuts = sorted(list(set(cuts)))
        for cut in cuts:
            dataCard = open('{0}/{1}_{2}_{3}.txt'.format(filePath, batch, variable, cut), 'w')
            cardMaker.card_producer(yields[cut], backgrounds, dataCard)
    else:
        dataCard = open('{0}/{1}_{2}.txt'.format(filePath, batch, variable), 'w')
        cardMaker.card_producer(yields, backgrounds, dataCard)
