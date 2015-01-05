#! /usr/bin/env python
import subprocess, shlex, datetime, copy, math, pickle
from multiprocessing import Process
import ROOT as r
import PlotProducer as pp

paramFile   = open('scripts/fcncParams.pkl', 'rb')
styles      = pickle.load(paramFile)
scales      = pickle.load(paramFile)
combos      = pickle.load(paramFile)

def FakeInfo():
    def __init__(self, cat, type, pct):
        self._category  = ''
        self._type      = ''
        self._highPtPct = pct 

if __name__ == '__main__':

    lepType = 'Electron' # 'Muon' or 'Electron'
    outDir  = 'plots/Fakes/{0}'.format(lepType)

    r.gROOT.SetBatch()
    r.gStyle.SetOptStat(0)

    canvas = r.TCanvas('canvas', 'canvas', 600, 450)

    legend = r.TLegend(0.6, 0.6, 0.89, 0.89)
    legend.SetFillColor(0)
    #legend.SetFillStyle(0)
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

    print ''.join(['& {0} '.format(cat) for cat in categories]),
    print '\\\\ \\hline'
    for fakeType in fakeTypes:
        print ' {0} '.format(fakeType),
        count = 0
        hCombined = 0.
        for category in categories:
            hist    = inFile.Get('{0}/{1}/h1_FakeablePt'.format(category, fakeType))

            if hist:
                pct     = hist.Integral(hist.FindBin(40), hist.GetNbinsX())/hist.Integral()
                fakeInfo.append(pct)
                print '& {0:.1%} '.format(pct),
            else:
                print '& -- ',
                continue

            pp.set_hist_style(hist, category, styles)
            legend.AddEntry(hist, category)

            if count == 0:
                hist.SetTitle('{0};p_{{T}}; Normalized Entries/Bin'.format(styles[fakeType][4]))
                hist.GetXaxis().SetRangeUser(0., 150.)
                hist.GetYaxis().SetTitleOffset(1.2)
                hist.SetMinimum(0.)
                hist.SetMaximum(1.25*hist.GetMaximum())

                hist.DrawNormalized()
                line.Draw('same')
            else:
                hist.DrawNormalized('same')

            if hCombined == 0.:
                hCombined = hist.Clone()
            else:
                hCombined.Add(hist)

            count += 1

        print '& {0:.1%} \\\\ \\hline'.format(hCombined.Integral(hCombined.FindBin(40.), hCombined.GetNbinsX())/hCombined.Integral())

        legend.Draw()
        canvas.Print('{0}/FakeablePt_{1}.png'.format(outDir, fakeType))
        legend.Clear()

        hCombined.DrawNormalized()
        canvas.Print('{0}/CombinedFakeablePt_{1}.png'.format(outDir, fakeType))

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
                        
                        nFinal = hYields.GetBinContent(9)
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
    dataFile    = r.TFile('data/fakeRates_DATA.root', 'OPEN')
    mcFile      = r.TFile('data/fakeRates_MC.root', 'OPEN')
    datasets    = ['ZJets', 'ttbar', 'QCD', 'WJets']

    for dataset in datasets:
        print dataset,
        h1_fakePt        = mcFile.Get('MC_truth_{0}/h1_{1}FakePt'.format(dataset, lepType))
        h1_fakePtLowJet  = mcFile.Get('MC_truth_{0}/h1_{1}FakePtLowJet'.format(dataset, lepType))
        h1_fakePtHighJet = mcFile.Get('MC_truth_{0}/h1_{1}FakePtHighJet'.format(dataset, lepType))

        if not h1_fakePtLowJet or not h1_fakePtHighJet: continue

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
    datasets    = ['QCD', 'WJets', 'ttbar', 'ZJets']

    mc_avg      = {}
    for variable in variables:
        binContent  = {}
        binError    = {}
        h1_avg      = r.TH1F()
        for i,dataset in enumerate(datasets):
            hist = mcFile.Get('MC_truth_{0}/h1_{1}{2}'.format(dataset, lepType, variable))
            if not hist: continue
            hist.SetTitle('#mu fake rates;{0};fake rate'.format(hist.GetXaxis().GetTitle()))
            pp.set_hist_style(hist, dataset, styles)

            binContent[dataset] = []
            binError[dataset]   = []
            for bin in range(hist.GetNbinsX()):
                binContent[dataset].append(hist.GetBinContent(bin))
                binError[dataset].append(hist.GetBinError(bin))

            if i == 0:
                hist.SetMinimum(0.)
                hist.SetMaximum(0.3)
                hist.Draw()

                h1_avg = hist.Clone()
            else:
                hist.Draw('same')


            legend.AddEntry(hist, styles[dataset][4])

        for bin in range(h1_avg.GetNbinsX()): 
            #h1_avg.SetBinContent(bin, binContent['WJets'][bin]) 
            #h1_avg.SetBinError(bin, binError['WJets'][bin]) 

            content, error = 0., 0.
            if variable == 'FakePtHighJet' or variable == 'FakePt': #post-selection
                if lepType == 'Muon':
                    error   = math.sqrt((0.63*binError['ttbar'][bin])**2 + (0.36*binError['WJets'][bin])**2 + (0.05*binError['QCD'][bin])**2) 
                    content = 0.63*binContent['ttbar'][bin] + 0.36*binContent['WJets'][bin] + 0.05*binContent['QCD'][bin] 
                if lepType == 'Electron':
                    error   = math.sqrt((0.26*binError['ttbar'][bin])**2 + (0.73*binError['WJets'][bin])**2 + (0.00*binError['QCD'][bin])**2) 
                    content = 0.26*binContent['ttbar'][bin] + 0.73*binContent['WJets'][bin] + 0.00*binContent['QCD'][bin] 
            #elif variable == 'FakePtLowJet': #pre-selection
            #    h1_avg.SetBinContent(bin, 0.021*binContent['ttbar'][bin] + 0.16*binContent['WJets'][bin] + 0.8*binContent['QCD'][bin]) 
            #    h1_avg.SetBinError(bin, math.sqrt((0.021*binError['ttbar'][bin])**2 + (0.16*binError['WJets'][bin])**2 + (0.8*binError['QCD'][bin])**2)) 
            else:
                #print binError['QCD'][bin], binError['ZJets'][bin]
                if binError['QCD'][bin] != 0 and binError['ZJets'][bin] != 0:
                    error   = 1./math.sqrt(1./binError['ZJets'][bin]**2 + 1./binError['QCD'][bin]**2)
                    content = (binContent['ZJets'][bin]/binError['ZJets'][bin]**2 + binContent['QCD'][bin]/binError['QCD'][bin]**2)*error**2
                else:
                    error   = 0.
                    content = 0.

            h1_avg.SetBinContent(bin, content)
            h1_avg.SetBinError(bin, error)

            #h1_avg.SetBinError(bin, 0.25*h1_avg.GetBinContent(bin+1))



        h1_avg.Draw('E2 SAME')
        h1_avg.SetLineColor(r.kGreen+3)
        h1_avg.SetMarkerColor(r.kGreen+3)
        h1_avg.SetMarkerStyle(20)
        h1_avg.SetFillColor(r.kGreen+3)
        h1_avg.SetFillStyle(3004)
        h1_avg.SetTitle('#mu fake rates (average);p_{{T}};fake rate')
        mc_avg[variable] = h1_avg.Clone()

        legend.AddEntry(h1_avg, 'Combined #pm 25%')
        legend.Draw()
        canvas.Print('{0}/{1}_overlays.png'.format(outDir, variable))
        legend.Clear()

    ### Overlay QCD from data and MC
    h1_qcd_MC    = mcFile.Get('MC_truth_QCD/h1_MuonFakePt')
    h1_qcd_Data  = dataFile.Get('QCD2l_DATA/h1_MuonFakePt')
    h1_qcd_MC.SetBit(r.TH1.kIsAverage)
    h1_qcd_Data.SetBit(r.TH1.kIsAverage)

    pp.set_hist_style(h1_qcd_MC, 'QCD', styles)
    pp.set_hist_style(h1_qcd_Data, 'DATA', styles)

    h1_qcd_MC.Draw('E')
    h1_qcd_Data.Draw('E SAME')

    legend.AddEntry(h1_qcd_MC, 'QCD (MC)')
    legend.AddEntry(h1_qcd_Data, 'QCD (Data)')
    legend.Draw()

    canvas.Print('{0}/{1}_overlays.png'.format(outDir, 'QCD_DataMC'))

    ### Overlay Z+jet from data and MC
    zJetsFile = r.TFile('data/fakeRates_ZJets.root', 'OPEN')

    #h1_zJets_MC      = zJetsFile.Get('ZPlusJet/h1_MuonFakePt')
    h1_zJets_MC      = mcFile.Get('ZPlusJet_ZJets/h1_MuonFakePt')
    h1_zJets_Data    = dataFile.Get('ZPlusJet_DATA/h1_MuonFakePt')
    h1_zJets_MC.SetBit(r.TH1.kIsAverage)
    h1_zJets_Data.SetBit(r.TH1.kIsAverage)

    pp.set_hist_style(h1_zJets_MC, 'ZJets', styles)
    pp.set_hist_style(h1_zJets_Data, 'DATA', styles)

    h1_zJets_MC.GetYaxis().SetRangeUser(0., 0.35)
    h1_zJets_MC.SetTitle(' #mu fake rates;p_{T};fake rate')

    h1_zJets_MC.Draw('E')
    h1_zJets_Data.Draw('E SAME')

    legend.Clear()
    legend.AddEntry(h1_zJets_MC, 'ZPlusJet (MC)')
    legend.AddEntry(h1_zJets_Data, 'ZPlusJet (Data)')
    legend.Draw()

    canvas.Print('{0}/{1}_overlays.png'.format(outDir, 'ZPlusJet_DataMC'))

    ### Overlay Z+jet and QCD from data
    h1_data_avg = h1_qcd_Data.Clone()
    h1_data_avg.Add(h1_zJets_Data)

    pp.set_hist_style(h1_qcd_Data, 'QCD', styles)
    pp.set_hist_style(h1_zJets_Data, 'ZJets', styles)
    pp.set_hist_style(mc_avg['FakePt'], 'AVG', styles)

    h1_qcd_Data.GetYaxis().SetRangeUser(0., 0.4)
    h1_qcd_Data.SetTitle(' #mu fake rates;p_{T};fake rate')

    h1_data_avg.GetYaxis().SetRangeUser(0., 0.4)
    h1_data_avg.GetXaxis().SetRangeUser(10., 39.)
    #h1_qcd_Data.Draw('E')
    #h1_zJets_Data.Draw('E SAME')
    #h1_data_avg.Draw('E')
    #mc_avg['FakePtHighJet'].Draw('E2')
    mc_avg['FakePt'].Draw('E2')
    mc_avg['FakePtLowJet'].Draw('E2 SAME')

    legend.Clear()
    #legend.AddEntry(h1_data_avg, 'Combined (Data)')
    #legend.AddEntry(mc_avg['FakePt'], 'Combined #geq 2 jets (MC)')
    #legend.AddEntry(mc_avg['FakePtLowJet'], 'Combined < 2 jets (MC)')
    legend.AddEntry(mc_avg['FakePt'], 'signal CR')
    legend.AddEntry(mc_avg['FakePtLowJet'], 'measurement CR')
    legend.Draw()

    canvas.Print('{0}/{1}_overlays.png'.format(outDir, 'DataVsMC_JetSplit'))

    ### Ratio between low and high jet mc fake rates
    h1_Ratio = mc_avg['FakePt'].Clone()
    h1_Ratio.Divide(mc_avg['FakePtLowJet'])
    #h1_data_avg.Add(h1_zJets_Data)

    pp.set_hist_style(h1_Ratio, 'AVG', styles)

    h1_Ratio.SetTitle(';p_{T};FR_{MC,signal}/FR_{MC,measured}')
    h1_Ratio.GetYaxis().SetRangeUser(0., 2.)
    h1_Ratio.GetYaxis().CenterTitle()
    h1_Ratio.GetYaxis().SetTitleSize(0.045)
    h1_Ratio.SetLineWidth(2)
    fitResult = h1_Ratio.Fit('pol0', 'S', '', 10, 100)

    h1_AvgBand = r.TH1F('h1_AvgBand', '', 1, 10., 100.)
    h1_AvgBand.SetLineWidth(0)
    h1_AvgBand.SetFillColor(r.kRed)
    h1_AvgBand.SetFillStyle(3004)
    h1_AvgBand.SetBinContent(1, fitResult.Parameters()[0])
    h1_AvgBand.SetBinError(1, 0.2*fitResult.Parameters()[0])

    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.1)
    r.gStyle.SetOptFit(1)
    #h1_data_avg.GetYaxis().SetRangeUser(0., 0.45)
    #h1_qcd_Data.Draw('E')
    #h1_zJets_Data.Draw('E SAME')
    h1_Ratio.Draw('E')
    h1_AvgBand.Draw('E2 SAME')

    canvas.Print('{0}/{1}.png'.format(outDir, 'Ratio'))


    ### Compare QCD jet multiplicity rates
    #h1_qcdJet_MC    = mcFile.Get('MC_truth_QCD/h1_MuonFakeJetMult')
    #h1_qcdJet_Data  = dataFile.Get('QCD2l/h1_MuonFakeJetMult')
    #h1_qcdJet_MC.SetBit(r.TH1.kIsAverage)
    #h1_qcdJet_Data.SetBit(r.TH1.kIsAverage)

    #pp.set_hist_style(h1_qcdJet_MC, 'QCD', styles)
    #pp.set_hist_style(h1_qcdJet_Data, 'DATA', styles)

    #h1_qcdJet_MC.Draw('E')
    #h1_qcdJet_Data.Draw('E SAME')

    #legend.AddEntry(h1_qcdJet_MC, 'QCD (MC)')
    #legend.AddEntry(h1_qcdJet_Data, 'QCD (Data)')
    #legend.Draw()

    #canvas.Print('{0}/{1}_overlays.png'.format(outDir, 'QCD_jets_DataMC'))

    ### Plot jet flavor for various backgrounds
    inFile = r.TFile('fakeEstimator/histos/20141130_234120.root')
    #for dataset in datasets:
