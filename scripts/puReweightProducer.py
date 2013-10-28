#! /usr/bin/env python
import ROOT as r
from decimal import *


r.gStyle.SetOptStat(0)
dataFile    = r.TFile('data/puHistograms_XS_71000.root', 'OPEN')
mcFile      = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20131024_004428.root', 'OPEN')
savePath    = 'plots'
canvas      = r.TCanvas()

h1_Data = dataFile.Get('pileup')
h1_MC   = mcFile.GetDirectory('inclusive/ZJets_M-50').Get('h1_SimVertexMultTrue')

r.TH1.SetDefaultSumw2(r.kTRUE)

h1_Data.Scale(1./h1_Data.Integral())
h1_MC.Scale(1./h1_MC.Integral())

h1_Data.SetLineColor(4)
h1_MC.SetLineColor(1)

h1_MC.SetMaximum(max(h1_Data.GetMaximum(), h1_MC.GetMaximum())*1.2)

h1_MC.SetAxisRange(0., 60.)
h1_Data.SetAxisRange(0., 60.)

#canvas.SetLogy()

h1_MC.Draw('HIST')
h1_Data.Draw('HIST SAME')

legend = r.TLegend(0.7, 0.75, 0.89, 0.89)
legend.SetFillStyle(4001)
legend.SetFillColor(0)
legend.SetTextSize(0.035)
legend.AddEntry(h1_Data, 'Run2012ABCD')
legend.AddEntry(h1_MC, 'S10 PU')
legend.Draw()

canvas.SaveAs(savePath+'/puPVMult.png')

outFile = r.TFile('data/puReweight_71000.root', 'RECREATE') 
canvas.SetGridy()

h1_PU  = r.TH1D("h1_PU", "2012ABCD PU reweight factors;PU;#omega_{PU}", 500, 0., 100.)
h1_PU.Divide(h1_Data, h1_MC)

for i in range(500): 
    h1_PU.SetBinError(i+1, 0.01)

h1_PU.SetMarkerStyle(33)
h1_PU.SetMarkerColor(r.kBlue)

line = r.TLine(0., 1., 100., 1.) 
line.SetLineWidth(2)
line.SetLineColor(r.kRed)

h1_PU.SetMaximum(4)
h1_PU.Draw('E0')
line.Draw('SAME')
canvas.SaveAs(savePath+'/puWeights_2012.png')

outFile.Write()
outFile.Close()
