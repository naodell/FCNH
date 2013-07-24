#! /usr/bin/env python
import subprocess
import ROOT as r

canvas  = r.TCanvas('canvas', 'canvas', 800, 600)
inFile  = r.TFile('fakeHistograms.root', 'OPEN')
outFile = r.TFile('fakeRates.root', 'RECREATE') 

elCategories = ['ele_v1', 'ele_v2', 'ele_v3', 'ele_v4'] 
muCategories = ['mu_v1']#, 'mu_v2'] 

fakeCategories = []
#fakeCategories.extend(elCategories)
fakeCategories.extend(muCategories)


canvas.SetGridx()
canvas.SetGridy()

for cat in fakeCategories:

    h1_NumerPt = inFile.GetDirectory(cat).Get("h1_NumerPt")
    h1_DenomPt = inFile.GetDirectory(cat).Get("h1_DenomPt")

    g_EffPt = r.TGraphAsymmErrors()
    g_EffPt.Divide(h1_NumerPt, h1_DenomPt)
    g_EffPt.SetName('g_FakePt_' + cat)
    g_EffPt.SetTitle(cat + ' fake rate;p_{T};#varepsilon')
    g_EffPt.Draw("AP")

    #canvas.Print("test.png")

    outFile.Add(g_EffPt)

    h1_NumerEta = inFile.GetDirectory(cat).Get("h1_NumerEta")
    h1_DenomEta = inFile.GetDirectory(cat).Get("h1_DenomEta")

    g_EffEta = r.TGraphAsymmErrors()
    g_EffEta.Divide(h1_NumerEta, h1_DenomEta)
    g_EffEta.SetName('g_FakeEta_' + cat)
    g_EffEta.SetTitle(cat + ' fake rate;#eta;#varepsilon')
    g_EffEta.Draw("AP")

    #canvas.Print("test.png")

    outFile.Add(g_EffEta)

outFile.Write()
outFile.Close()

