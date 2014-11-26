#! /usr/bin/env python
import subprocess, shlex, datetime, copy, math, pickle
from multiprocessing import Process
import ROOT as r

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

    now         = datetime.datetime.now()
    currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

    r.gROOT.SetBatch()

    ### One component of the systematic uncertainty is from extrapolating to
    ### fake rates above 40 GeV.  To assess how important this is, first calculate
    ### the percentage of events with a contribution from fakes with pt > 40 GeV
    inFile = r.TFile('fcncAnalysis/combined_histos/fcnh_cut4_2012_20141120_144325.root', 'OPEN')
    categories  = ['ss_mumu', 'ss_ee', 'ss_emu', '3l_inclusive']
    fakeTypes   = ['muFakes', 'eFakes', 'llFakes']
    fakeInfo    = []
    lumi        = 19700

    print ''.join(['| *{0}* '.format(cat) for cat in categories]),
    print '|'
    for fakeType in fakeTypes:
        print '| {0} '.format(fakeType),
        for category in categories:
            hist    = inFile.Get('{0}/{1}/h1_FakeablePt'.format(category, fakeType))

            if hist:
                pct     = hist.Integral(hist.FindBin(40), hist.GetNbinsX())/hist.Integral()
                fakeInfo.append(pct)
                print '| {0:.1%} '.format(pct),
            else:
                print '| -- ',

        print '|'

    print '\n\n'

    ### Additionally we would like to know the contribution from different
    ### processes in the signal region.  It seems only ttbar, W+jet, Z+jets, and, maybe,
    ### QCD will be important.
    datasets    = ['WJets', 'ZJets', 'ttbarHad', 'ttbarLep', 'QCD', 'BG', 'DATA']
    sampleSets  = {
                'DATA':['DATA_MUON', 'DATA_ELECTRON', 'DATA_MUEG'],
                'WJets':['W1JetsToLNu', 'W2JetsToLNu', 'W3JetsToLNu', 'W4JetsToLNu'],
                'ZJets': ['ZJets_M-50', 'ZJets_M-10To50'],
                'ttbar':['ttbarHad', 'ttabarLep'],
                #'ttbarHad':['ttbarHad'],
                #'ttbarLep':['ttbarLep'],
                'QCD':['QCD_20_MU', 'QCD_20-30_EM', 'QCD_30-80_EM', 'QCD_80-170_EM', 'QCD_170-250_EM', 'QCD_250-350_EM', 'QCD_350_EM']
              }
        
    totals = {}
    for category in categories:
        totals[category] = {}

        totals[category]['BG'] = 0
        for dataset,samples in sampleSets.iteritems():
            subtotal = 0
            for sample in samples:
                hYields = inFile.Get('{0}/{1}/h1_YieldByCut'.format(category, sample))

                if hYields:

                    if dataset != 'DATA':
                        nInit   = hYields.GetBinContent(1)
                        hYields.Scale(lumi*scales['2012'][sample]/nInit)
                    
                    nFinal = hYields.GetBinContent(9)
                    subtotal += nFinal

            totals[category][dataset]   = subtotal
            if dataset != 'DATA':
                totals[category]['BG']      += subtotal


    print '| ** ',
    for dataset in datasets:
        print '| *{0}* '.format(dataset),
    print '|'
        
    for category in categories:
        print '| {0} '.format(category),
        for dataset in datasets:
            if dataset not in ['DATA', 'BG']:
                if totals[category]['BG'] == 0: 
                    print '| -- ',
                else:
                    print '| {0:.2%} '.format(totals[category][dataset]/totals[category]['BG']),
            else:
                print '| {0:.2f} '.format(totals[category][dataset]),
        print '|'
            
    ### Check dependency of fake rate on jet multiplicity
    lepType     = 'Muon' # 'Muon' or 'Electron'
    datasets    = ['ZJets', 'ttbar', 'QCD', 'WJets'g]
    rateFile    = r.TFile('data/fakeRates_TEST.root')
    outDir      = 'plots/Fakes/{0}/'.format(now)

    for dataset in datasets:
        h1_fakePt = rateFile.Get('MC_Truth_{0}/h1_MuonFakePt'.format())
        h1_fakePtLowJet  = rateFile.Get('MC_Truth_{0}/h1_MuonFakePtLowJet'.format())
        h1_fakePtHighJet = rateFile.Get('MC_Truth_{0}/h1_MuonFakePtHighJet'.format())


    exit()
