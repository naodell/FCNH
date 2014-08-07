#! /usr/bin/env python
import subprocess, shlex, time, pickle, math, sys
from array import array
import ROOT as r

seed = int(str(time.time())[6:10])

'''
Carries out a random grid search to optimize rectangular cuts in 2D plane.
'''

class yields():
    def __init(self, category):
        self._category  = ''
        self._samples   = []
        self._yields    = []

#def card_producer(yields):
    # 

if len(sys.argv) > 0:
    batch   = sys.argv[1]
else:
    print 'A batch and ratio type must be specified.  Otherwise, do some hacking so this thing knows about your inputs.'
    exit()

categories  = ['ee', 'emu', 'mumu']
datasets    = ['fcnh', 'irr', 'Fakes', 'QFlips']

dataDict = {}
dataDict['Fakes']     = ['muFakes', 'eFakes', 'llFakes'] # Fakes
dataDict['QFlips']    = ['QFlips'] # electron charge misID
dataDict['irr']       = ['WZJets3LNu', 'ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau', 'ttZ', 'ttW', 'ttG']) # Irreducible backgrounds
dataDict['fcnh']      = ['FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_ZZ_t', 'FCNC_ZZ_tbar', 'FCNC_TauTau_t', 'FCNC_TauTau_tbar'] # signal

varList = ['metVsHt']

hFile = r.TFile('fcncAnalysis/combined_histos/fcnh_cut3_2012_{0}.root'.format(batch), 'OPEN')

canvas = r.TCanvas('canvas', 'canvas', 800, 600)
canvas.SetGridx()
canvas.SetGridy()

# Scale factors
paramFile = open('scripts/fcncParams.pkl', 'rb')
scales    = pickle.load(paramFile)

for category in categories:

    # Get histograms for different samples 
    histDict = {}
    for dataset in datasets:
        histDict[dataset] = []

        sumBG = r.TH2D()
        for i,sample in enumerate(dataDict[dataset]):
            hist = hFile.GetDirectory(category + '/' + sample).Get('h2_' + var)

            if not hist:
                continue

            #print sample
            hist.Scale(scales['2012'][sample])

            if i is 0:
                sumBG = hist
            else:
                sumBG.Add(hist)

    xLow    = sumBG.GetXaxis().GetXmin() 
    xHigh   = sumBG.GetXaxis().GetXmax() 
    yLow    = sumBG.GetYaxis().GetXmin() 
    yHigh   = sumBG.GetYaxis().GetXmax() 

    xDiv    = sumBG.GetNbinsX()
    yDiv    = sumBG.GetNbinsY()

    #print 'Using random seed {}'.format(seed)
    rndNum = r.TRandom3(seed)
    effSig  = []
    effBG   = []

    maxSig      = [0, 0, 0]
    maxSigCuts  = [0, 0]

    maxEffRatio = [0, 0, 0]
    maxEffCuts  = [0, 0]

    xCutHigh = xDiv
    yCutHigh = yDiv
    count = 0
    #for xCutLow in range(1,xDiv):
    for xCutLow in range(1,xDiv):
        for yCutLow in range(1,yDiv):

            numSig  = sumSig.Integral(xCutLow, xCutHigh, yCutLow, yCutHigh)
            denSig  = sumSig.Integral()
            numBG   = sumBG.Integral(xCutLow, xCutHigh, yCutLow, yCutHigh)
            denBG   = sumBG.Integral()

            effSig.append(numSig/denSig) 
            effBG.append(numBG/denBG) 

            #print (xCutLow, xCutHigh), (yCutLow, yCutHigh), (numSig, denSig), (numBG, denBG), significance, numSig*denBG/(denSig*numBG)
            if numSig != 0 or numBG != 0:
                significance = numSig/math.sqrt(numSig + numBG)
                if significance > maxSig[0]:
                    maxSig      = [significance, effSig[count], effBG[count]]
                    maxSigCuts  = [xHigh*(xCutLow/float(xDiv)), xHigh*(xCutHigh/float(xDiv)), yHigh*(yCutLow/float(yDiv)), yHigh*(yCutHigh/float(yDiv))]

                effRatio = numSig/denSig - numBG/denBG                
                print effRatio, xCutLow, yCutLow
                if effRatio > maxEffRatio[0]:
                    maxEffRatio = [effRatio, effSig[count], effBG[count]]
                    maxEffCuts  = [xHigh*(xCutLow/float(xDiv)), xHigh*(xCutHigh/float(xDiv)), yHigh*(yCutLow/float(yDiv)), yHigh*(yCutHigh/float(yDiv))]

            count += 1

    print maxSig, maxSigCuts
    print maxEffRatio, maxEffCuts

    sumSig.Draw('colz')
    canvas.SaveAs('plots/sumSig_' + var + '_test.png')
    sumBG.Draw('colz')
    canvas.SaveAs('plots/sumBG_' + var + '_test.png')

    line  = r.TLine(0., 0., 1.1, 1.1)
    line.SetLineColor(r.kRed)
    line.SetLineWidth(2)

    g_eff = r.TGraph(len(effSig), array('d', effBG), array('d', effSig))
    g_eff.SetTitle(var + ' ROC;#varepsilon_{BG};#varepsilon_{signal}')
    g_eff.SetMarkerStyle(21)
    g_eff.SetMarkerColor(r.kBlue)
    g_eff.Draw('AP')
    line.Draw('same')

    canvas.SaveAs('plots/eff_' + var + '_test.png')

