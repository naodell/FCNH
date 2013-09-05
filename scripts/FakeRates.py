#! /usr/bin/env python
from array import array
import subprocess
import ROOT as r

canvas  = r.TCanvas('canvas', 'canvas', 800, 600)
inFile  = r.TFile('histos/fakeHistograms.root', 'OPEN')
outFile = r.TFile('histos/fakeRates.root', 'RECREATE') 
path    = '~/work/plots/fakes/rates'

fakeCategories = []
fakeCategories.extend(['ele_v1', 'ele_v2', 'ele_v3', 'ele_v4'])
fakeCategories.extend(['mu_v1', 'mu_v2']) 

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

h2_Fakes2D  = []
ptBins      = [10., 15., 20., 25., 30., 35., 50., 80., 120., 250.];
etaBins     = [-2.5, -2., -1.479, -1., 0., 1., 1.479, 2., 2.5];

for i,cat in enumerate(fakeCategories):

    print cat
    canvas.SetGridx()
    canvas.SetGridy()

    # pt fakes
    h1_NumerPt = inFile.GetDirectory(cat + '/DATA').Get('h1_NumerPt')
    h1_DenomPt = inFile.GetDirectory(cat + '/DATA').Get('h1_DenomPt')

    g_EffPt = r.TGraphAsymmErrors()
    g_EffPt.Divide(h1_NumerPt, h1_DenomPt)
    g_EffPt.SetName('g_FakePt_' + cat)
    g_EffPt.SetTitle(cat + ' fake rate;p_{T};#varepsilon')
    g_EffPt.Draw('AP')

    canvas.Print(path + '/' + cat + '_pt.png')

    outFile.Add(g_EffPt)

    # eta fakes
    h1_NumerEta = inFile.GetDirectory(cat + '/DATA').Get('h1_NumerEta')
    h1_DenomEta = inFile.GetDirectory(cat + '/DATA').Get('h1_DenomEta')

    g_EffEta = r.TGraphAsymmErrors()
    g_EffEta.Divide(h1_NumerEta, h1_DenomEta)
    g_EffEta.SetName('g_FakeEta_' + cat)
    g_EffEta.SetTitle(cat + ' fake rate;#eta;#varepsilon')
    g_EffEta.Draw('AP')

    canvas.Print(path + '/' + cat + '_eta.png')

    outFile.Add(g_EffEta)

    # 2D (pt and eta) fakes
    canvas.SetLogy()
    canvas.SetGridx(0)
    canvas.SetGridy(0)

    h2_Numer = inFile.GetDirectory(cat + '/DATA').Get('h2_NumerPtVsEta')
    h2_Denom = inFile.GetDirectory(cat + '/DATA').Get('h2_DenomPtVsEta')

    h2_Fakes2D.append(r.TH2D('h2_Fakes_' + cat, 'fake rates ' + cat +';#eta;p_{T}', 5, array('f', etaBins), 6, array('f', ptBins)))
    #h2_Fakes2D[i].Divide(h2_Numer, h2_Denom)
    #h2_Fakes2D[i].Draw('colz text')

    #canvas.Print(path + '/' + cat + '_2D.png')
    canvas.SetLogy(0)


outFile.Write()
outFile.Close()
