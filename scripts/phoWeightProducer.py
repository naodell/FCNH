#! /usr/bin/env python
import ROOT as r
from PlotProducer import *

#path = '~/afs/public_html/higgs/2011/November/3/weights/gamma'
path = '~/afs/public_html/higgs/TEST/weights/gamma'
doBGSubtraction = True
doPtWeights     = True 
doPVWeights     = True 
doJetWeights    = True 
doTests         = False

 ### Input files ###
#muonFile     = r.TFile('../HiggsAnalyzer/histos/archive/Nov1/muon_NoPVWeights_fullPhoton.root', 'OPEN')

gammaFile    = r.TFile('../HiggsAnalyzer/histos/higgsHistograms_Photons.root', 'OPEN')
muonFile     = r.TFile('../HiggsAnalyzer/histos/higgsHistograms_Muons.root', 'OPEN')
electronFile = r.TFile('../HiggsAnalyzer/histos/higgsHistograms_Electrons.root', 'OPEN')

r.gROOT.SetStyle('Plain')
r.gStyle.SetOptStat(0)
canvas    = r.TCanvas()
canvas.SetLogy()
canvas.SetGridy()
canvas.SetGridx()

histPtGamma    = gammaFile.GetDirectory('Lepton').Get('h1_diLeptonQt_Photon')
histPtMuon     = muonFile.GetDirectory('Lepton').Get('h1_diLeptonQt_DATA')
histPtElectron = electronFile.GetDirectory('Lepton').Get('h1_diLeptonQt_DATA')

#histPtGamma    = gammaFile.GetDirectory('Weights').Get('h1_QtForWeights_Photon')
#histPtMuon     = muonFile.GetDirectory('Weights').Get('h1_QtForWeights_DATA')
#histPtElectron = electronFile.GetDirectory('Weights').Get('h1_QtForWeights_DATA')

histPVMuGamma  = muonFile.GetDirectory('Misc').Get('h1_pvMult_PhotonJets')
histPVMuon     = muonFile.GetDirectory('Misc').Get('h1_pvMult_DATA')
histPVEGamma   = electronFile.GetDirectory('Misc').Get('h1_pvMult_PhotonJets')
histPVElectron = electronFile.GetDirectory('Misc').Get('h1_pvMult_DATA')

histJetMuGamma  = muonFile.GetDirectory('Jet').Get('h1_jetMult_PhotonJets')
histJetMuon     = muonFile.GetDirectory('Jet').Get('h1_jetMult_DATA')
histJetEGamma   = electronFile.GetDirectory('Jet').Get('h1_jetMult_PhotonJets')
histJetElectron = electronFile.GetDirectory('Jet').Get('h1_jetMult_DATA')

 ### background subtraction ###

if doBGSubtraction:
    LUMI_MU = 1.0514*(.2151 + .9302 + .3709 + .6630)
    LUMI_E  = 1.053*(.2151 + .7982 + .3132 + .6332)
    bgList   = ['ttbar', 'tW', 'WW', 'WZ', 'ZZ']

    muPlotter  = PlotProducer('../HiggsAnalyzer/histos/higgsHistograms_Muons.root', LUMI_MU)
    #muHistList = muPlotter.ScaleHistsByDict(muPlotter.GetHistList('Weights', 'QtForWeights', bgList))
    muHistList = muPlotter.ScaleHistsByDict(muPlotter.GetHistList('Lepton', 'diLeptonQt', bgList))
    muSumHist  = muPlotter.GetSumHist(muHistList)
    muSumHist.SetTitle('muon BG')
    muSumHist.Draw('E3')
    muPlotter.MakeSavePath(path)
    canvas.SaveAs(path+'/MuonBG.png')
    histPtMuon.Add(muSumHist, -1)

    ePlotter  = PlotProducer('../HiggsAnalyzer/histos/higgsHistograms_Electrons.root', LUMI_E)
    eHistList = ePlotter.ScaleHistsByDict(ePlotter.GetHistList('Lepton', 'diLeptonQt', bgList))
    eSumHist  = ePlotter.GetSumHist(eHistList)
    eSumHist.SetTitle('electron BG')
    eSumHist.Draw('E3')
    canvas.SaveAs(path+'/ElectronBG.png')
    histPtMuon.Add(eSumHist, -1)

 ### sanity check ###
histPtMuon.SetLineColor(r.kBlue)
histPtElectron.SetLineColor(r.kRed)
histPtMuon.SetMaximum(max(histPtMuon.GetMaximum(), histPtGamma.GetMaximum())*5)
histPtMuon.Draw('HIST')
histPtElectron.Draw('HIST SAME')
histPtGamma.Draw('HIST SAME')

legend = r.TLegend(0.70,0.60,0.87,0.86)
legend.SetFillColor(0)
legend.SetLineWidth(0)
legend.AddEntry(histPtGamma, '#gamma + jet')
legend.AddEntry(histPtMuon, '#mu#mu')
legend.AddEntry(histPtElectron, 'ee')
legend.Draw()

canvas.SaveAs(path+'/gammaJetQt.png')

outFile = r.TFile('../HiggsAnalyzer/data/gammaReweight.root', 'RECREATE') 

 ### mass spectrum ###

muGammaMass         = r.TH1D(muonFile.GetDirectory('Lepton').Get('h1_diLeptonMass_DATA'))
muGammaMass.SetName('h1_diMuonMass')

eGammaMass         = r.TH1D(electronFile.GetDirectory('Lepton').Get('h1_diLeptonMass_DATA'))
eGammaMass.SetName('h1_diElectronMass')

 ### pt weights ###
if doPtWeights:
    muGammaPtReweightHist = histPtMuon.Clone()
    muGammaPtReweightHist.SetName('h1_muGammaPtWeight')
    muGammaPtReweightHist.Divide(histPtMuon, histPtGamma)
    muGammaPtReweightHist.Draw()
    canvas.SaveAs(path+'/muGammaPtWeights.png')

    eGammaPtReweightHist = histPtElectron.Clone()
    eGammaPtReweightHist.SetName('h1_eGammaPtWeight')
    eGammaPtReweightHist.Divide(histPtElectron, histPtGamma)
    eGammaPtReweightHist.Draw()
    canvas.SaveAs(path+'/eGammaPtWeights.png')


 ### PV weights ### 
if doPVWeights:
    canvas.SetLogy(0)
    muGammaPVReweightHist = r.TH1D("h1_muGammaPVWeight", "#gamma PV reweight factors (muons)", 25, 0.5, 25.5)
    muGammaPVReweightHist.Divide(histPVMuon, histPVMuGamma)
    muGammaPVReweightHist.Draw()
    canvas.SaveAs(path+'/muGammaPVWeights.png')

    eGammaPVReweightHist = r.TH1D("h1_eGammaPVWeight", "#gamma PV reweight factors (electrons)", 25, 0.5, 25.5)
    eGammaPVReweightHist.Divide(histPVElectron, histPVEGamma)
    eGammaPVReweightHist.Draw()
    canvas.SaveAs(path+'/eGammaPVWeights.png')

 ### jet weights ###
if doJetWeights:
    bin1 = histJetMuon.GetBinContent(1)/histJetMuGamma.GetBinContent(1)
    bin2 = histJetMuon.GetBinContent(2)/histJetMuGamma.GetBinContent(2)
    print 'muon ~ 0 jet: %f\t 1 jet: %f' % (bin1, bin2)

    bin1 = histJetElectron.GetBinContent(1)/histJetEGamma.GetBinContent(1)
    bin2 = histJetElectron.GetBinContent(2)/histJetEGamma.GetBinContent(2)
    print 'muon ~ 0 jet: %f\t 1 jet: %f' % (bin1, bin2)

### TESTS!!! ###
if doTests:
    muGammaPtReweightHist = testPtMuon.Clone()
    muGammaPtReweightHist.Divide(testPtMuon, testPtGamma)
    muGammaPtReweightHist.SetAxisRange(0., 270.)
    muGammaPtReweightHist.SetLineWidth(2)
    muGammaPtReweightHist.Draw()
    canvas.SaveAs(path+'/testPtWeights.png')

outFile.Write()
outFile.Close()
