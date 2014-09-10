#! /usr/bin/env python
import subprocess, shlex, datetime, copy, math
from multiprocessing import Process


import ROOT as r

if __name__ == '__main__':

    r.gROOT.SetBatch()

    now         = datetime.datetime.now()
    currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

    lepType     = 'Muon' # 'Muon' or 'Electron'
    datasets    = ['DATA', 'ZJets', 'ttbar', 'QCD', 'WJets', 'DATA_QCD2l']
    dataFiles   = {}
    for dataset in datasets:
        dataFiles[dataset] = r.TFile('data/fakeRates_{0}.root'.format(dataset), 'OPEN')

    h1_fakeRates = {}

    ### Rates measured in data ###
    h1_fakeRates['ZJets_DATA']  = dataFiles['DATA'].Get('ZPlusJet/h1_{0}FakePt'.format(lepType))
    h1_fakeRates['ZJets_DATA'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['ZJets_DATA'].SetLineColor(r.kBlue)
    h1_fakeRates['ZJets_DATA'].SetLineWidth(2)

    h1_fakeRates['QCD2l_DATA']  = dataFiles['DATA_QCD2l'].Get('QCD2l/h1_{0}FakePt'.format(lepType))
    h1_fakeRates['QCD2l_DATA'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['QCD2l_DATA'].SetLineColor(r.kRed)
    h1_fakeRates['QCD2l_DATA'].SetLineWidth(2)

    h1_fakeRates['Combined_DATA'] = h1_fakeRates['ZJets_DATA'].Clone()
    h1_fakeRates['Combined_DATA'].SetTitle('fake rates measured in data;p_{T};fake rate')
    h1_fakeRates['Combined_DATA'].GetYaxis().SetTitleOffset(1.5)
    h1_fakeRates['Combined_DATA'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['Combined_DATA'].Add(h1_fakeRates['QCD2l_DATA'])
    h1_fakeRates['Combined_DATA'].SetLineColor(r.kGreen)
    h1_fakeRates['Combined_DATA'].SetLineWidth(0)
    h1_fakeRates['Combined_DATA'].SetFillColor(r.kGreen)
    h1_fakeRates['Combined_DATA'].SetFillStyle(3002)
    h1_fakeRates['Combined_DATA'].SetMarkerColor(r.kGreen)
    h1_fakeRates['Combined_DATA'].SetMarkerStyle(22)
    h1_fakeRates['Combined_DATA'].GetXaxis().SetRange(0, 6)
    h1_fakeRates['Combined_DATA'].GetYaxis().SetRangeUser(0., 0.5)

    for bin in range(h1_fakeRates['Combined_DATA'].GetNbinsX()):
        h1_fakeRates['Combined_DATA'].SetBinError(bin, math.sqrt(pow(h1_fakeRates['Combined_DATA'].GetBinError(bin), 2) + pow(0.25*h1_fakeRates['Combined_DATA'].GetBinContent(bin), 2)))

    canvas = r.TCanvas('canvas', 'canvas', 650, 700)
    canvas.SetLeftMargin(0.15)

    h1_fakeRates['Combined_DATA'].Draw('E2')
    h1_fakeRates['QCD2l_DATA'].Draw('SAME')
    h1_fakeRates['ZJets_DATA'].Draw('SAME')

    legend = r.TLegend(0.5,0.65,.89,0.89)
    legend.SetFillColor(0)
    legend.SetLineWidth(1)
    legend.SetLineColor(1)
    legend.SetTextSize(0.03)

    legend.AddEntry(h1_fakeRates['Combined_DATA'], 'Combined #pm 25% syst')
    legend.AddEntry(h1_fakeRates['QCD2l_DATA'], 'b#bar{b}')
    legend.AddEntry(h1_fakeRates['ZJets_DATA'], 'Z+jet')
    legend.Draw()

    canvas.SaveAs('plots/Fakes/fake_syst_DATA.pdf')

    ### Z+jets rates ###
    h1_fakeRates['ZJets_ZJets'] = dataFiles['ZJets'].Get('ZPlusJet/h1_{0}FakePt'.format(lepType))
    h1_fakeRates['ZJets_ZJets'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['ZJets_ZJets'].SetLineColor(r.kRed)
    h1_fakeRates['ZJets_ZJets'].SetLineWidth(2)

    h1_fakeRates['MC_truth_ZJets']  = dataFiles['ZJets'].Get('MC_truth/h1_{0}FakePt'.format(lepType))
    h1_fakeRates['MC_truth_ZJets'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['MC_truth_ZJets'].SetLineColor(r.kGreen+4)
    h1_fakeRates['MC_truth_ZJets'].SetLineWidth(2)

    h1_fakeRates['Combined_ZJets'] = h1_fakeRates['ZJets_DATA'].Clone()
    h1_fakeRates['Combined_ZJets'].SetTitle('fake rates from Z+jets;p_{T};fake rate')
    h1_fakeRates['Combined_ZJets'].GetYaxis().SetTitleOffset(1.5)
    h1_fakeRates['Combined_ZJets'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['Combined_ZJets'].Add(h1_fakeRates['ZJets_ZJets'])
    h1_fakeRates['Combined_ZJets'].Add(h1_fakeRates['MC_truth_ZJets'])
    h1_fakeRates['Combined_ZJets'].SetLineColor(r.kGreen)
    h1_fakeRates['Combined_ZJets'].SetLineWidth(0)
    h1_fakeRates['Combined_ZJets'].SetFillColor(r.kGreen)
    h1_fakeRates['Combined_ZJets'].SetFillStyle(3002)
    h1_fakeRates['Combined_ZJets'].SetMarkerColor(r.kGreen)
    h1_fakeRates['Combined_ZJets'].SetMarkerStyle(22)
    h1_fakeRates['Combined_ZJets'].GetXaxis().SetRange(0, 6)
    h1_fakeRates['Combined_ZJets'].GetYaxis().SetRangeUser(0., 0.5)

    for bin in range(h1_fakeRates['Combined_ZJets'].GetNbinsX()):
        h1_fakeRates['Combined_ZJets'].SetBinError(bin, math.sqrt(pow(h1_fakeRates['Combined_ZJets'].GetBinError(bin), 2) + pow(0.25*h1_fakeRates['Combined_ZJets'].GetBinContent(bin), 2)))

    h1_fakeRates['Combined_ZJets'].Draw('E2')
    h1_fakeRates['ZJets_DATA'].Draw('SAME')
    h1_fakeRates['ZJets_ZJets'].Draw('SAME')
    h1_fakeRates['MC_truth_ZJets'].Draw('SAME')


    legend.Clear()
    legend.AddEntry(h1_fakeRates['Combined_ZJets'], 'Combined #pm 25%')
    legend.AddEntry(h1_fakeRates['ZJets_DATA'], 'Z+jet (Data)')
    legend.AddEntry(h1_fakeRates['ZJets_ZJets'], 'Z+jet (MC)')
    legend.AddEntry(h1_fakeRates['MC_truth_ZJets'], 'Z+jet (MC-truth)')
    legend.Draw()

    canvas.SaveAs('plots/Fakes/fake_syst_ZJets.pdf')

    ### MC rates ###
    h1_fakeRates['MC_truth_QCD']    = dataFiles['QCD'].Get('MC_truth/h1_{0}FakePt'.format(lepType))
    h1_fakeRates['MC_truth_QCD'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['MC_truth_QCD'].SetLineColor(r.kRed-4)
    h1_fakeRates['MC_truth_QCD'].SetLineWidth(2)

    h1_fakeRates['MC_truth_WJets']  = dataFiles['WJets'].Get('MC_truth/h1_{0}FakePt'.format(lepType))
    h1_fakeRates['MC_truth_WJets'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['MC_truth_WJets'].SetLineColor(r.kViolet)
    h1_fakeRates['MC_truth_WJets'].SetLineWidth(2)

    h1_fakeRates['MC_truth_ttbar']  = dataFiles['ttbar'].Get('MC_truth/h1_{0}FakePt'.format(lepType))
    h1_fakeRates['MC_truth_ttbar'].SetBit(r.TH1.kIsAverage)
    h1_fakeRates['MC_truth_ttbar'].SetLineColor(r.kBlue)
    h1_fakeRates['MC_truth_ttbar'].SetLineWidth(2)

    #h1_fakeRates['Combined_MC'] = h1_fakeRates['MC_truth_ZJets'].Clone()
    #h1_fakeRates['Combined_MC'].Add(h1_fakeRates['MC_truth_WJets'])

    #h1_fakeRates['Combined_MC'] = h1_fakeRates['MC_truth_WJets'].Clone()
    #h1_fakeRates['Combined_MC'].SetBit(r.TH1.kIsAverage)
    #h1_fakeRates['Combined_MC'].Add(h1_fakeRates['MC_truth_QCD'])
    #h1_fakeRates['Combined_MC'].Add(h1_fakeRates['MC_truth_ttbar'])

    h1_fakeRates['Combined_MC'] = h1_fakeRates['MC_truth_WJets'].Clone()
    h1_fakeRates['Combined_MC'].SetBit(r.TH1.kIsAverage)

    for bin in range(1,7):
        nWJets = h1_fakeRates['MC_truth_WJets'].GetBinContent(bin)
        nttbar = h1_fakeRates['MC_truth_ttbar'].GetBinContent(bin)
        nQCD   = h1_fakeRates['MC_truth_QCD'].GetBinContent(bin)
        h1_fakeRates['Combined_MC'].SetBinContent(bin, 0.05*nWJets + 0.25*nttbar + 0.7*nQCD)
        h1_fakeRates['Combined_MC'].SetBinError(bin, 0.25*h1_fakeRates['Combined_MC'].GetBinContent(bin))

    h1_fakeRates['Combined_MC'].SetTitle('fake rates from MC truth;p_{T};fake rate')
    h1_fakeRates['Combined_MC'].SetLineColor(r.kGreen)
    h1_fakeRates['Combined_MC'].SetLineWidth(0)
    h1_fakeRates['Combined_MC'].SetFillColor(r.kGreen)
    h1_fakeRates['Combined_MC'].SetFillStyle(3002)
    h1_fakeRates['Combined_MC'].SetMarkerColor(r.kGreen)
    h1_fakeRates['Combined_MC'].SetMarkerStyle(22)
    h1_fakeRates['Combined_MC'].GetXaxis().SetRange(0, 6)
    h1_fakeRates['Combined_MC'].GetYaxis().SetRangeUser(0., 0.5)
    h1_fakeRates['Combined_MC'].GetYaxis().SetTitleOffset(1.5)

    h1_fakeRates['Combined_MC'].Draw('E2')
    h1_fakeRates['MC_truth_QCD'].Draw('SAME')
    #h1_fakeRates['MC_truth_ZJets'].Draw('SAME')
    h1_fakeRates['MC_truth_WJets'].Draw('SAME')
    h1_fakeRates['MC_truth_ttbar'].Draw('SAME')


    legend.Clear()
    legend.AddEntry(h1_fakeRates['Combined_MC'], 'Combined #pm 25%')
    legend.AddEntry(h1_fakeRates['MC_truth_QCD'], 'QCD')
    #legend.AddEntry(h1_fakeRates['MC_truth_ZJets'], 'Z+jet')
    legend.AddEntry(h1_fakeRates['MC_truth_WJets'], 'W+jet')
    legend.AddEntry(h1_fakeRates['MC_truth_ttbar'], 't#bar{t}')
    legend.Draw()

    canvas.SaveAs('plots/Fakes/fake_syst_MC.pdf')

    ### Data vs. MC ###

    h1_fakeRates['Combined_MC'].SetLineColor(r.kRed)
    h1_fakeRates['Combined_MC'].SetFillColor(r.kRed)
    h1_fakeRates['Combined_MC'].SetFillStyle(3005)
    h1_fakeRates['Combined_MC'].SetMarkerColor(r.kRed)

    h1_fakeRates['Combined_DATA'].SetLineColor(r.kBlue)
    h1_fakeRates['Combined_DATA'].SetFillColor(r.kBlue)
    h1_fakeRates['Combined_DATA'].SetFillStyle(3004)
    h1_fakeRates['Combined_DATA'].SetMarkerColor(r.kBlue)

    h1_fakeRates['Combined_MC'].SetTitle('Data Vs. MC fake rate;p_{T};fake rate')

    
    bError = [0., 0.]
    for bin in range(h1_fakeRates['Combined_MC'].GetNbinsX()):
        bError[0] = h1_fakeRates['Combined_DATA'].GetBinError(bin)
        bError[1] = h1_fakeRates['Combined_MC'].GetBinError(bin)

        h1_fakeRates['Combined_DATA'].SetBinError(bin, 0.01)
        h1_fakeRates['Combined_MC'].SetBinError(bin, 0.01)

    h1_fakeRates['Combined_MIX'] = h1_fakeRates['Combined_MC'].Clone()
    h1_fakeRates['Combined_MIX'].Add(h1_fakeRates['Combined_DATA'])
    h1_fakeRates['Combined_MIX'].SetLineColor(r.kGreen)
    h1_fakeRates['Combined_MIX'].SetFillColor(r.kGreen)
    h1_fakeRates['Combined_MIX'].SetFillStyle(3002)
    h1_fakeRates['Combined_MIX'].SetMarkerStyle(0)
    h1_fakeRates['Combined_MIX'].SetMarkerColor(r.kGreen)

    for bin in range(h1_fakeRates['Combined_MC'].GetNbinsX()):
        h1_fakeRates['Combined_DATA'].SetBinError(bin, 0.25*h1_fakeRates['Combined_DATA'].GetBinContent(bin))
        h1_fakeRates['Combined_MC'].SetBinError(bin, 0.25*h1_fakeRates['Combined_MC'].GetBinContent(bin))


    for bin in range(h1_fakeRates['Combined_MIX'].GetNbinsX()):
        h1_fakeRates['Combined_MIX'].SetBinError(bin, abs(h1_fakeRates['Combined_MC'].GetBinContent(bin) - h1_fakeRates['Combined_DATA'].GetBinContent(bin))/2.)

    h1_fakeRates['Combined_MC'].Draw('E2')
    h1_fakeRates['Combined_DATA'].Draw('E2 SAME')
    h1_fakeRates['Combined_MIX'].Draw('E2 SAME')

    legend.Clear()
    legend.AddEntry(h1_fakeRates['Combined_MC'], 'Combined MC #pm 25%')
    legend.AddEntry(h1_fakeRates['Combined_DATA'], 'Combined Data #pm 25%')
    legend.AddEntry(h1_fakeRates['Combined_MIX'], 'Data - MC')
    legend.Draw()

    canvas.SaveAs('plots/Fakes/fake_syst_DataVsMC.pdf')
