#! /usr/bin/env python
import subprocess, shlex, datetime, copy
from multiprocessing import Process
import ROOT as r

now         = datetime.datetime.now()
currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

iso     = [4, 5, 6, 7, 8, 9]
batches = ['20140630_160143', '20140630_155842', '20140630_152157', '20140630_152028', '20140630_151832', '20140630_145739']
isoDict = dict(zip(batches, iso))
hists   = {}

files = []
for batch in batches:
    file = r.TFile('fcncAnalysis/combined_histos/fcnh_cut3_2012_{0}.root'.format(batch), 'OPEN')
    files.append(file)

    h1_JetMult_mu = file.GetDirectory('ss_mumu/muFakes').Get('h1_AllJetMult')
    h1_JetMult_ll = file.GetDirectory('ss_mumu/llFakes').Get('h1_AllJetMult')
    h1_JetMult_comb = h1_JetMult_mu.Clone()
    h1_JetMult_comb.Add(h1_JetMult_ll)

    h1_JetMult_comb.SetLineColor(r.kBlue+(isoDict[batch]-4))
    h1_JetMult_comb.SetLineWidth(2)

    hists[isoDict[batch]] = []
    hists[isoDict[batch]].append(h1_JetMult_mu)
    hists[isoDict[batch]].append(h1_JetMult_ll)
    hists[isoDict[batch]].append(h1_JetMult_comb)

h1_JetMult_DATA = files[0].GetDirectory('ss_mumu/DATA_MUON').Get('h1_AllJetMult')

canvas = r.TCanvas('canvas', 'canvas', 650, 700)

h1_JetMult_DATA.SetMaximum(1.3*h1_JetMult_DATA.GetMaximum())
h1_JetMult_DATA.Draw('E1')
#hists[4][2].Draw('HIST SAME')
hists[5][2].Draw('HIST SAME')
hists[6][2].Draw('HIST SAME')
hists[7][2].Draw('HIST SAME')
hists[8][2].Draw('HIST SAME')
#hists[9][2].Draw('HIST SAME')

canvas.SaveAs('plots/TEST.png')

