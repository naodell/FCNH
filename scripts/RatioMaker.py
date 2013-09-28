#! /usr/bin/env python
from array import array
from PlotProducer import *
import subprocess
import ROOT as r

canvas  = r.TCanvas('canvas', 'canvas', 800, 600)
inFile  = r.TFile('histos/fcnh_cut1_2012_20130822_183318.root', 'OPEN')
inPath  = 'inclusive/ttbar'
#inFile  = r.TFile('histos/fcncHistograms_cut1.root', 'OPEN')
#inPath  = 'inclusive/TEST'
outFile = r.TFile('histos/bEff_ttbar_2012.root', 'RECREATE') 
outPath = 'plots/bTag_eff'

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

ratioDict   = {'bTagEff':('BTruthNumerPt', 'BTruthDenomPt'), 'bMistagEff':('BMistagNumerPt', 'BMistagDenomPt')}

canvas.SetGridx()
canvas.SetGridy()

for key,value in ratioDict.iteritems():

    print 'Producing ratio plot {}'.format(key)

    h1_Numer = inFile.GetDirectory(inPath).Get('h1_{0}'.format(value[0]))
    h1_Denom = inFile.GetDirectory(inPath).Get('h1_{0}'.format(value[1]))

    g_Eff = r.TGraphAsymmErrors()
    g_Eff.Divide(h1_Numer, h1_Denom)
    g_Eff.SetName('g_{0}'.format(key))
    g_Eff.SetTitle('{0};p_{{T}};#varepsilon'.format(key))
    #g_Eff.Draw('ACP')

    #canvas.Print('{0}/{1}.png'.format(outPath, key))

    outFile.Add(g_Eff)

outFile.Write()
outFile.Close()
