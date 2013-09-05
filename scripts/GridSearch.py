#! /usr/bin/env python
import subprocess, shlex, time, pickle, math
from array import array
import ROOT as r

seed = int(str(time.time())[6:10])

'''
Carries out a random grid search to optimize rectangular cuts in 2D plane.
'''

nToys   = 5000

bgList  = []
bgList.extend(['ZJets', 'ZJets_M-10To50', 'WJets']) # V+jets
bgList.extend(['WWJets2L2Nu', 'ZZJets2L2Nu', 'ZZJets2L2Q', 'WZJets2L2Q']) # Diboson to 2l + X
bgList.extend(['WZJets3LNu']) # WZ to 3l+nu
bgList.extend(['ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau']) # ZZ to 4l
bgList.extend(['tW', 'ttbar', 'ttZ', 'ttW']) # Top

sigList = ['FCNC_M125']

varList = ['metVsHt', 'TrileptonMVsDileptonMOS']
catList = ['3l_inclusive']

hFile = r.TFile('histos/fcncHistograms_cut1.root', 'OPEN')

canvas = r.TCanvas('canvas', 'canvas', 800, 600)
canvas.SetGridx()
canvas.SetGridy()

# Scale factors
paramFile = open('scripts/fcncParams.pkl', 'rb')
scales    = pickle.load(paramFile)

for var in varList:

    # Get signal histograms
    sumSig = r.TH2D()
    for i,sample in enumerate(sigList):
        hist = hFile.GetDirectory(catList[0] + '/' + sample).Get('h2_' + var)

        if not hist:
            continue

        hist.Scale(scales['2012'][sample])

        if i is 0:
            sumSig = hist
        else:
            sumSig.Add(hist)

    # Get background sums
    sumBG = r.TH2D()
    for i,sample in enumerate(bgList):
        hist = hFile.GetDirectory(catList[0] + '/' + sample).Get('h2_' + var)

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

    for i in range(nToys):

        rArrayX = sorted([rndNum.Rndm(), rndNum.Rndm()])
        rArrayY = sorted([rndNum.Rndm(), rndNum.Rndm()])

        xCutLow    = int(math.floor(xDiv*rArrayX[0])) 
        #xCutHigh   = xDiv
        xCutHigh   = int(math.floor(xDiv*rArrayX[1]))
        yCutLow    = int(math.floor(yDiv*rArrayY[0]))
        #yCutHigh   = yDiv
        yCutHigh   = int(math.floor(yDiv*rArrayY[1]))

        effSig.append(sumSig.Integral(xCutLow, xCutHigh, yCutLow, yCutHigh)/sumSig.Integral()) 
        effBG.append(sumBG.Integral(xCutLow, xCutHigh, yCutLow, yCutHigh)/sumBG.Integral()) 

        if effSig[i] != 0 or effBG[i] != 0:
            significance = effSig[i]/math.sqrt(effSig[i] + effBG[i])

            if significance > maxSig[0]:
                maxSig      = [significance, effSig[i], effBG[i]]
                maxSigCuts  = [xHigh*(xCutLow/float(xDiv)), xHigh*(xCutHigh/float(xDiv)), yHigh*(yCutLow/float(yDiv)), yHigh*(yCutHigh/float(yDiv))]

                print maxSig, maxSigCuts

    #print maxSig, maxSigCuts

    g_eff = r.TGraph(len(effSig), array('d', effBG), array('d', effSig))
    g_eff.SetTitle(var + ' ROC;#varepsilon_{BG};#varepsilon_{signal}')
    g_eff.SetMarkerStyle(21)
    g_eff.SetMarkerColor(r.kBlue)
    g_eff.Draw('AP')

    canvas.SaveAs('../plots/TEST/eff_' + var + 'test.png')

