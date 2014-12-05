#! /usr/bin/env python
import subprocess, shlex, datetime, copy, math, pickle
from multiprocessing import Process
import ROOT as r
import PlotProducer as pp

paramFile   = open('scripts/fcncParams.pkl', 'rb')
scales      = pickle.load(paramFile)
styles      = pickle.load(paramFile)
combos      = pickle.load(paramFile)

def FakeInfo():
    def __init__(self, cat, type, pct):
        self._category  = ''
        self._type      = ''
        self._highPtPct = pct 

if __name__ == '__main__':

    outDir      = 'plots/Fakes/TEST'

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)

    canvas = r.TCanvas('canvas', 'canvas', 600, 450)

    legend = r.TLegend(0.65, 0.6, 0.89, 0.89)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    ### One component of the systematic uncertainty is from extrapolating to
    ### fake rates above 40 GeV.  To assess how important this is, first calculate
    ### the percentage of events with a contribution from fakes with pt > 40 GeV
    inFile = r.TFile('fcncAnalysis/combined_histos/fcnh_cut4_2012_20141120_144325.root', 'OPEN')
    categories  = ['ss_mumu', 'ss_ee', 'ss_emu', '3l_inclusive']
    fakeTypes   = ['muFakes', 'eFakes']#, 'llFakes']
    fakeInfo    = []
    lumi        = 19700

    line = r.TLine(40., 40., 0., 1.)
    line.SetLineWidth(2)

    print ''.join(['| *{0}* '.format(cat) for cat in categories]),
    print '|'
    for fakeType in fakeTypes:
        print '| {0} '.format(fakeType),
        count = 0
        for category in categories:
            hist    = inFile.Get('{0}/{1}/h1_FakeablePt'.format(category, fakeType))

            if hist:
                pct     = hist.Integral(hist.FindBin(40), hist.GetNbinsX())/hist.Integral()
                fakeInfo.append(pct)
                print '| {0:.1%} '.format(pct),
            else:
                print '| -- ',
                continue

            pp.set_hist_style(hist, category, styles)
            legend.AddEntry(hist, category)

            if count == 0:
                hist.SetTitle('{0};p_{{T}}; Normalized Entries/Bin'.format(styles[fakeType][4]))
                hist.GetXaxis().SetRangeUser(0., 150.)
                hist.SetMinimum(0.)
                hist.SetMaximum(1.25*hist.GetMaximum())

                hist.DrawNormalized()
                line.Draw('same')
            else:
                hist.DrawNormalized('same')

            count += 1

        legend.Draw()
        canvas.Print('{0}/FakeablePt_{1}.png'.format(outDir, fakeType))
        print '|'
        legend.Clear()

    print '\n\n'

    ### Additionally we would like to know the contribution from different
    ### processes in the signal region.  It seems only ttbar, W+jet, Z+jets, and, maybe,
    ### QCD will be important.
    datasets    = ['WJets', 'ZJets', 'ttbar', 'QCD']#, 'DATA']
    sampleSets  = {
                'DATA':['DATA_MUON', 'DATA_ELECTRON', 'DATA_MUEG'],
                'WJets':['W1JetsToLNu', 'W2JetsToLNu', 'W3JetsToLNu', 'W4JetsToLNu'],
                'ZJets': ['ZJets_M-50', 'ZJets_M-10To50'],
                'ttbar':['ttbarHad', 'ttabarLep'],
                'QCD':['QCD_20_MU', 'QCD_20-30_EM', 'QCD_30-80_EM', 'QCD_80-170_EM', 'QCD_170-250_EM', 'QCD_250-350_EM', 'QCD_350_EM']
              }
        
    totals = {}
    for category in categories:
        totals[category] = {}

        totals[category]['BG'] = 0
        for dataset in datasets:
            subtotal = 0
            if dataset in sampleSets.keys():
                for sample in sampleSets[dataset]: 
                    hYields = inFile.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))
                    if hYields:

                        if dataset != 'DATA':
                            nInit   = hYields.GetBinContent(1)
                            hYields.Scale(lumi*scales['2012'][sample]/nInit)
                        
                        nFinal = hYields.GetBinContent(6)
                        subtotal += nFinal
            else:
                hYields = inFile.Get('{0}/{1}/h1_YieldByCut'.format(category, dataset))
                if hYields:

                    if dataset[:4] != 'DATA':
                        nInit   = hYields.GetBinContent(1)
                        hYields.Scale(lumi*scales['2012'][dataset]/nInit)
                    
                    nFinal = hYields.GetBinContent(6)
                    subtotal += nFinal

            totals[category][dataset] = subtotal

            if dataset != 'DATA': 
                totals[category]['BG'] += subtotal

    datasets.append('BG')

    print ' ',
    for dataset in datasets:
        print '& {0} '.format(dataset),
    print '\\\\ \hline'
        
    for category in categories:
        print ' {0} '.format(category),
        for dataset in datasets:
            if dataset not in ['DATA', 'BG']:
                if totals[category]['BG'] == 0: 
                    print '& -- ',
                else:
                    print '& {0:.1%} '.format(totals[category][dataset]/totals[category]['BG']),
            else:
                print '& {0:.2f} '.format(totals[category][dataset]),
        print '\\\\'
    print '\hline'
            
    ### Check dependency of fake rate on jet multiplicity
    lepType     = 'Muon' # 'Muon' or 'Electron'
    mcFile    = r.TFile('data/fakeRates_TEST.root', 'OPEN')
    datasets    = ['ZJets', 'ttbar', 'QCD', 'WJets']

    for dataset in datasets:
        h1_fakePt        = mcFile.Get('MC_truth_{0}/h1_{1}FakePt'.format(dataset, lepType))
        h1_fakePtLowJet  = mcFile.Get('MC_truth_{0}/h1_{1}FakePtLowJet'.format(dataset, lepType))
        h1_fakePtHighJet = mcFile.Get('MC_truth_{0}/h1_{1}FakePtHighJet'.format(dataset, lepType))

        h1_fakePt.SetTitle('#mu fake rates;p_{T};fake rate')

        h1_fakePt.SetLineColor(r.kBlack)
        h1_fakePt.SetLineWidth(2)
        h1_fakePtLowJet.SetLineColor(r.kBlue)

        h1_fakePtLowJet.SetLineWidth(2)
        h1_fakePtHighJet.SetLineColor(r.kRed)
        h1_fakePtHighJet.SetLineWidth(2)

        maxVal = max([h1_fakePt.GetMaximum(), h1_fakePtLowJet.GetMaximum(), h1_fakePtHighJet.GetMaximum()])
        h1_fakePt.SetMaximum(1.25*maxVal)

        h1_fakePt.Draw()
        h1_fakePtLowJet.Draw('same')
        h1_fakePtHighJet.Draw('same')

        legend.AddEntry(h1_fakePt, 'inclusive')
        legend.AddEntry(h1_fakePtLowJet, 'N_{jets} < 2')
        legend.AddEntry(h1_fakePtHighJet, 'N_{jets} >= 2')
        legend.Draw()

        canvas.Print('{0}/FakePt_{1}.png'.format(outDir, dataset))
        legend.Clear()

    ### Overlay QCD, W+jets, and ttbar
    variables   = ['FakePt', 'FakePtLowJet', 'FakePtHighJet']#, 'FakeJetMult']
    datasets    = ['QCD', 'WJets', 'ttbar']
    for variable in variables:
        binContent  = {}
        h1_avg      = r.TH1F()
        for i,dataset in enumerate(datasets):
            hist = mcFile.Get('MC_truth_{0}/h1_{1}{2}'.format(dataset, lepType, variable))
            hist.SetTitle('#mu fake rates;{0};fake rate'.format(hist.GetXaxis().GetTitle()))
            pp.set_hist_style(hist, dataset, styles)

            binContent[dataset] = []
            for bin in range(hist.GetNbinsX()):
                binContent[dataset].append(hist.GetBinContent(bin))

            if i == 0:
                hist.SetMinimum(0.)
                hist.SetMaximum(0.3)
                hist.Draw()

                h1_avg = hist.Clone()
                h1_avg.SetTitle('#mu fake rates (average);p_{{T}};fake rate')
            else:
                hist.Draw('same')

            legend.AddEntry(hist, styles[dataset][4])

        for bin in range(h1_avg.GetNbinsX()): 
            if variable == 'FakePtHighJet':
                h1_avg.SetBinContent(bin, 0.63*binContent['ttbar'][bin] + 0.36*binContent['WJets'][bin] + 0.1*binContent['QCD'][bin]) #post-selection
            else:
                h1_avg.SetBinContent(bin, 0.021*binContent['ttbar'][bin] + 0.16*binContent['WJets'][bin] + 0.8*binContent['QCD'][bin]) #pre-selection
            h1_avg.SetBinError(bin+1, 0.25*h1_avg.GetBinContent(bin+1))

        h1_avg.SetLineColor(r.kGreen)
        h1_avg.SetFillColor(r.kGreen)
        h1_avg.SetFillStyle(3004)
        h1_avg.Draw('E2 SAME')
        legend.AddEntry(h1_avg, 'Combined #pm 25%')
        
        legend.Draw()
        canvas.Print('{0}/{1}_overlays.png'.format(outDir, variable))

        legend.Clear()

    ### Overlay QCD from data and MC
    dataFile = r.TFile('data/fakeRates.root', 'OPEN')

    h1_qcd_MC   = mcFile.Get('MC_truth_QCD/h1_MuonFakePt')
    g_qcd_Data = dataFile.Get('QCD2l/g_MuonFakePt')

    pp.set_hist_style(h1_qcd_MC, 'QCD', styles)
    pp.set_hist_style(g_qcd_Data, 'DATA', styles)

    h1_qcd_MC.Draw()
    g_qcd_Data.Draw('P SAME')

    legend.AddEntry(h1_qcd_MC, 'QCD (MC)')
    legend.AddEntry(g_qcd_Data, 'QCD (Data)')
    legend.Draw()

    canvas.Print('{0}/{1}_overlays.png'.format(outDir, 'QCD_DataMC'))

    ### Overlay Z+jet from data and MC
    zJetsFile = r.TFile('data/fakeRates_ZJets.root', 'OPEN')

    g_zJets_MC      = zJetsFile.Get('ZPlusJet/g_MuonFakePt')
    g_zJets_Data    = dataFile.Get('ZPlusJet/g_MuonFakePt')

    pp.set_hist_style(g_zJets_MC, 'ZJets', styles)
    pp.set_hist_style(g_zJets_Data, 'DATA', styles)

    g_zJets_MC.GetYaxis().SetRangeUser(0., 0.35)
    g_zJets_MC.SetTitle(' #mu fake rates;p_{T};fake rate')

    g_zJets_MC.Draw('AP')
    g_zJets_Data.Draw('P SAME')

    legend.Clear()
    legend.AddEntry(g_zJets_MC, 'ZPlusJet (MC)')
    legend.AddEntry(g_zJets_Data, 'ZPlusJet (Data)')
    legend.Draw()

    canvas.Print('{0}/{1}_overlays.png'.format(outDir, 'ZPlusJet_DataMC'))

    ### Overlay Z+jet and QCD from data
    pp.set_hist_style(g_qcd_Data, 'QCD', styles)
    pp.set_hist_style(g_zJets_Data, 'ZJets', styles)

    g_qcd_Data.GetYaxis().SetRangeUser(0., 0.35)
    g_qcd_Data.SetTitle(' #mu fake rates;p_{T};fake rate')

    g_qcd_Data.Draw('AP')
    g_zJets_Data.Draw('P SAME')

    legend.Clear()
    legend.AddEntry(g_qcd_Data, 'QCD (Data)')
    legend.AddEntry(g_zJets_Data, 'ZPlusJet (Data)')
    legend.Draw()

    canvas.Print('{0}/{1}_overlays.png'.format(outDir, 'ZJetsVsQCD_DATA'))

    ### Plot jet flavor for various backgrounds
    inFile = r.TFile('fakeEstimator/histos/20141130_234120.root')
    #for dataset in datasets:
