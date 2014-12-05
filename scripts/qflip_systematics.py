#! /usr/bin/env python
import subprocess, shlex, datetime, copy, math
from multiprocessing import Process


import ROOT as r

if __name__ == '__main__':

    r.gROOT.SetBatch()

    now         = datetime.datetime.now()
    currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

    datasets = ['DATA', 'ZJets']#, 'ttbar']
    dataFiles = {}
    for dataset in datasets:
        dataFiles[dataset] = r.TFile('data/electronQMisID_{0}.root'.format(dataset), 'OPEN')

    g_QFlips = {}

    ### Rates measured in data ###
    g_QFlips['DATA_BB']  = dataFiles['DATA'].Get('inclusive/g_DielectronMisQ_BB')
    g_QFlips['DATA_BB'].SetLineColor(r.kBlue)
    g_QFlips['DATA_BB'].SetFillColor(0)
    g_QFlips['DATA_BB'].SetLineWidth(2)

    g_QFlips['DATA_BE']  = dataFiles['DATA'].Get('inclusive/g_DielectronMisQ_BE')
    g_QFlips['DATA_BE'].SetLineColor(r.kBlue)
    g_QFlips['DATA_BE'].SetFillColor(0)
    g_QFlips['DATA_BE'].SetLineWidth(2)

    g_QFlips['DATA_EE']  = dataFiles['DATA'].Get('inclusive/g_DielectronMisQ_EE')
    g_QFlips['DATA_EE'].SetLineColor(r.kBlue)
    g_QFlips['DATA_EE'].SetFillColor(0)
    g_QFlips['DATA_EE'].SetLineWidth(2)

    ### Rates measured in Z+jets MC ###
    g_QFlips['ZJets_BB']  = dataFiles['ZJets'].Get('inclusive/g_DielectronMisQ_BB')
    g_QFlips['ZJets_BB'].SetLineColor(r.kRed)
    g_QFlips['ZJets_BB'].SetFillColor(0)
    g_QFlips['ZJets_BB'].SetLineWidth(2)

    g_QFlips['ZJets_BE']  = dataFiles['ZJets'].Get('inclusive/g_DielectronMisQ_BE')
    g_QFlips['ZJets_BE'].SetLineColor(r.kRed)
    g_QFlips['ZJets_BE'].SetFillColor(0)
    g_QFlips['ZJets_BE'].SetLineWidth(2)

    g_QFlips['ZJets_EE']  = dataFiles['ZJets'].Get('inclusive/g_DielectronMisQ_EE')
    g_QFlips['ZJets_EE'].SetLineColor(r.kRed)
    g_QFlips['ZJets_EE'].SetFillColor(0)
    g_QFlips['ZJets_EE'].SetLineWidth(2)

    ### Rates measured using MC truth Z+jets ###
    #test  = dataFiles['ZJets'].GetDirectory('inclusive')
    #test.ls()
    #exit()

    g_QFlips['ZJets_MCInZ_BB']  = dataFiles['ZJets'].Get('inclusive/g_EleQMisIDInZ_MC_1')
    g_QFlips['ZJets_MCInZ_BB'].SetLineColor(r.kViolet)
    g_QFlips['ZJets_MCInZ_BB'].SetFillColor(0)
    g_QFlips['ZJets_MCInZ_BB'].SetLineWidth(2)

    g_QFlips['ZJets_MCInZ_BE']  = dataFiles['ZJets'].Get('inclusive/g_EleQMisIDInZ_MC_2')
    g_QFlips['ZJets_MCInZ_BE'].SetLineColor(r.kViolet)
    g_QFlips['ZJets_MCInZ_BE'].SetFillColor(0)
    g_QFlips['ZJets_MCInZ_BE'].SetLineWidth(2)

    g_QFlips['ZJets_MCInZ_EE']  = dataFiles['ZJets'].Get('inclusive/g_EleQMisIDInZ_MC_3')
    g_QFlips['ZJets_MCInZ_EE'].SetLineColor(r.kViolet)
    g_QFlips['ZJets_MCInZ_EE'].SetFillColor(0)
    g_QFlips['ZJets_MCInZ_EE'].SetLineWidth(2)

    g_QFlips['ZJets_MCOutZ_BB']  = dataFiles['ZJets'].Get('inclusive/g_EleQMisIDOutZ_MC_1')
    g_QFlips['ZJets_MCOutZ_BB'].SetLineColor(r.kGreen)
    g_QFlips['ZJets_MCOutZ_BB'].SetFillColor(0)
    g_QFlips['ZJets_MCOutZ_BB'].SetLineWidth(2)

    g_QFlips['ZJets_MCOutZ_BE']  = dataFiles['ZJets'].Get('inclusive/g_EleQMisIDOutZ_MC_2')
    g_QFlips['ZJets_MCOutZ_BE'].SetLineColor(r.kGreen)
    g_QFlips['ZJets_MCOutZ_BE'].SetFillColor(0)
    g_QFlips['ZJets_MCOutZ_BE'].SetLineWidth(2)

    g_QFlips['ZJets_MCOutZ_EE']  = dataFiles['ZJets'].Get('inclusive/g_EleQMisIDOutZ_MC_3')
    g_QFlips['ZJets_MCOutZ_EE'].SetLineColor(r.kGreen)
    g_QFlips['ZJets_MCOutZ_EE'].SetFillColor(0)
    g_QFlips['ZJets_MCOutZ_EE'].SetLineWidth(2)

    canvas = r.TCanvas('canvas', 'canvas', 650, 700)
    canvas.SetLeftMargin(0.15)

    legend = r.TLegend(0.55,0.21,.89,0.39)
    legend.SetFillColor(0)
    legend.SetLineWidth(1)
    legend.SetLineColor(1)
    legend.SetTextSize(0.03)

    g_QFlips['ZJets_MCInZ_BB'].SetMaximum(0.0008)
    g_QFlips['ZJets_MCInZ_BB'].SetMinimum(0.)
    g_QFlips['ZJets_MCInZ_BB'].SetTitle('Charge flips (|#eta| < 0.8);p_{T};charge flip rate')
    g_QFlips['ZJets_MCInZ_BB'].GetYaxis().SetTitleOffset(1.2)

    g_QFlips['ZJets_MCInZ_BB'].Draw('AL')
    g_QFlips['ZJets_BB'].Draw('L SAME')
    g_QFlips['DATA_BB'].Draw('L SAME')

    legend.AddEntry(g_QFlips['DATA_BB'], 'Data')
    legend.AddEntry(g_QFlips['ZJets_BB'], 'Z+jets')
    legend.AddEntry(g_QFlips['ZJets_MCInZ_BB'], 'Z+jets (MC truth)')
    legend.Draw()

    canvas.SaveAs('plots/qflip_syst_DataMC_BB.pdf')

    g_QFlips['ZJets_MCInZ_BE'].SetMaximum(0.006)
    g_QFlips['ZJets_MCInZ_BE'].SetMinimum(0.)
    g_QFlips['ZJets_MCInZ_BE'].SetTitle('Charge flips (0.8 < |#eta| < 1.479);p_{T};charge flip rate')
    g_QFlips['ZJets_MCInZ_BE'].GetYaxis().SetTitleOffset(1.2)

    g_QFlips['ZJets_MCInZ_BE'].Draw('AL')
    g_QFlips['ZJets_BE'].Draw('L SAME')
    g_QFlips['DATA_BE'].Draw('L SAME')
    legend.Draw()

    canvas.SaveAs('plots/qflip_syst_DataMC_BE.pdf')

    g_QFlips['ZJets_MCInZ_EE'].SetMaximum(0.012)
    g_QFlips['ZJets_MCInZ_EE'].SetMinimum(0.)
    g_QFlips['ZJets_MCInZ_EE'].SetTitle('Charge flips (1.479 < |#eta| < 2.1);p_{T};charge flip rate')
    g_QFlips['ZJets_MCInZ_EE'].GetYaxis().SetTitleOffset(1.2)

    g_QFlips['ZJets_MCInZ_EE'].Draw('AL')
    g_QFlips['ZJets_EE'].Draw('L SAME')
    g_QFlips['DATA_EE'].Draw('L SAME')
    legend.Draw()

    canvas.SaveAs('plots/qflip_syst_DataMC_EE.pdf')

    g_QFlips['ZJets_MCInZ_BB'].Draw('AL')
    g_QFlips['ZJets_MCOutZ_BB'].Draw('L SAME')

    legend.Clear()
    legend.AddEntry(g_QFlips['ZJets_MCInZ_BB'], 'Z+jets MC on-shell')
    legend.AddEntry(g_QFlips['ZJets_MCOutZ_BB'], 'Z+jets MC off-shell')
    legend.Draw()

    canvas.SaveAs('plots/qflip_syst_MCtruth_BB.pdf')

    g_QFlips['ZJets_MCInZ_BE'].Draw('AL')
    g_QFlips['ZJets_MCOutZ_BE'].Draw('L SAME')
    legend.Draw()

    canvas.SaveAs('plots/qflip_syst_MCtruth_BE.pdf')

    g_QFlips['ZJets_MCInZ_EE'].Draw('AL')
    g_QFlips['ZJets_MCOutZ_EE'].Draw('L SAME')
    legend.Draw()

    canvas.SaveAs('plots/qflip_syst_MCtruth_EE.pdf')
