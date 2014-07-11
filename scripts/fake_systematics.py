#! /usr/bin/env python
import subprocess, shlex, datetime, copy
from multiprocessing import Process
import ROOT as r

now         = datetime.datetime.now()
currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

iso     = [4, 5, 6, 7, 8, 9]
batches = ['20140709_011007', '20140709_010159', '20140709_005529', '20140709_005047', '20140709_004655', '20140709_004516']
isoDict = dict(zip(batches, iso))

files = {}
for batch in batches:
    file = r.TFile('fcncAnalysis/combined_histos/fcnh_cut1_2012_{0}.root'.format(batch), 'OPEN')
    files[batch] = file

varNames = ['Met', 'Lepton1Pt', 'Lepton2Pt', 'AllJetMult']
for varName in varNames:
    hists   = {}
    for batch in batches:
        h_mu = files[batch].GetDirectory('ss_mumu/muFakes').Get('h1_' + varName)
        h_ll = files[batch].GetDirectory('ss_mumu/llFakes').Get('h1_' + varName)
        h_comb = h_mu.Clone()
        h_comb.Add(h_ll)

        h_comb.SetLineColor(r.kBlue+(isoDict[batch]-4))
        h_comb.SetLineWidth(2)

        hists[isoDict[batch]] = []
        hists[isoDict[batch]].append(h_mu)
        hists[isoDict[batch]].append(h_ll)
        hists[isoDict[batch]].append(h_comb)

    h_DATA = files['20140709_011007'].GetDirectory('ss_mumu/DATA_MUON').Get('h1_' + varName)

    canvas = r.TCanvas('canvas', 'canvas', 650, 700)

    h_DATA.SetMaximum(1.3*h_DATA.GetMaximum())
    h_DATA.Draw('E1')
    #hists[4][2].Draw('HIST SAME')
    hists[5][2].Draw('HIST SAME')
    hists[6][2].Draw('HIST SAME')
    hists[7][2].Draw('HIST SAME')
    hists[8][2].Draw('HIST SAME')
    #hists[9][2].Draw('HIST SAME')

    canvas.SaveAs('plots/IsoStudy_{0}.png'.format(varName))


