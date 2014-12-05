#! /usr/bin/env python
import subprocess, shlex, datetime, copy
from multiprocessing import Process
import ROOT as r

r.gStyle.SetOptStat(0)

now         = datetime.datetime.now()
currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

histFile_pre    = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_20141202_021714.root', 'OPEN')
histFile_post   = r.TFile('fcncAnalysis/combined_histos/fcnh_cut3_2012_20141202_021714.root', 'OPEN')

canvas = r.TCanvas('canvas', 'canvas', 650, 700)

legend = r.TLegend(0.55, 0.6, 0.89, 0.89)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.03)

histStyles = [(r.kBlue, '0.2 < Iso_{rel} < 0.5'), 
                (r.kGreen+3, '0.5 < Iso_{rel} < 0.8'), 
                (r.kRed, '0.8 < Iso_{rel} < 1.1'), 
                (r.kMagenta, '1.1 < Iso_{rel} < 1.4')] 

histNames = ['FakeableIsoRelVsPt', 'FakeableIsoRelVsMET', 'FakeableIsoRelVsJetMultiplicity']
for histName in histNames:
    hist = histFile_pre.Get('ss_mumu/muFakes/h2_{0}'.format(histName))
    
    projections = []
    for i,bin in enumerate(range(3,hist.GetNbinsY(),3)):
        projection = hist.ProjectionX('h1_{0}_{1}'.format(histName, bin), bin, bin+2)
        projections.append(projection)

        projection.SetLineWidth(2)
        projection.SetLineColor(histStyles[i][0])

        if i == 0:
            projection.DrawNormalized()
        else:
            projection.DrawNormalized('same')

        legend.AddEntry(projection, histStyles[i][1])

    legend.Draw()
    canvas.SaveAs('plots/IsoStudy_{0}_test.png'.format(histName))
    legend.Clear()


