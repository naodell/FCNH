#! /usr/bin/env python
from array import array
from PlotProducer import *
import subprocess
import ROOT as r

canvas  = r.TCanvas('canvas', 'canvas', 800, 600)
inFile  = r.TFile('histos/fcnh_cut1_2012_20130821_164730.root', 'OPEN')
inPath  = 'inclusive/ttbar'
outFile = r.TFile('histos/bEff_ttbar_2012.root', 'RECREATE') 
outPath = '~/work/plots/bTag_eff'

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

h2_Fakes2D  = []
ptBins      = [10., 15., 20., 25., 30., 35., 50., 80., 120., 250.];
etaBins     = [-2.5, -2., -1.479, -1., 0., 1., 1.479, 2., 2.5];

ratioDict   = {'bTagEff':('BTruthNumerPt', 'BTruthDenomPt'), 'bMistagEff':('BMistagNumerPt', 'BTruthDenomPt')}

canvas.SetGridx()
canvas.SetGridy()

for key,value in ratioDict.iteritems():

    print 'Producing ratio plot {}'.format(key)
    # pt fakes
    h1_Numer = inFile.GetDirectory(inPath).Get('h1_{0}'.format(value[0]))
    h1_Denom = inFile.GetDirectory(inPath).Get('h1_{0}'.format(value[1]))

    g_Eff = r.TGraphAsymmErrors()
    g_Eff.Divide(h1_Numer, h1_Denom)
    g_Eff.SetName('g_FakePt_' + cat)
    g_Eff.SetTitle(cat + ' fake rate;p_{T};#varepsilon')
    g_Eff.Draw('AP')

    canvas.Print('{0}/{1}.png'.format(outPath, key)

    outFile.Add(g_Eff)

outFile.Write()
outFile.Close()
