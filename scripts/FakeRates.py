#! /usr/bin/env python
from RatioMaker import *

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(r.kTRUE)
r.TH2.SetDefaultSumw2(r.kTRUE)

if len(sys.argv) > 1:
    batch       = sys.argv[1]
    selection   = sys.argv[2]
else:
    print 'You forgot to specify what file to run over and whether you want DATA or MC.'
    exit()

doCombination = False

inFile  = 'fakeEstimator/histos/{0}.root'.format(batch)
outFile = 'data/fakeRates_TEST.root'

fakeCategories = []
datasets = {}
bgType = ''
if selection == 'DATA':
    datasets['QCD2l']       = ['DATA']
    datasets['ZPlusJet']    = ['DATA']
    bgType = 'PROMPT'

    fakeCategories.append('QCD2l')
    fakeCategories.append('ZPlusJet')
elif selection == 'MC':
    datasets['QCD2l']       = ['QCD', 'ttbar', 'ZJets', 'WJets']
    datasets['ZPlusJet']    = ['ttbar', 'ZJets']
    datasets['SameSign']    = ['QCD', 'ttbar', 'ZJets', 'WJets']
    datasets['MC_truth']    = ['QCD', 'ttbar', 'ZJets', 'WJets']

    fakeCategories.append('QCD2l')
    fakeCategories.append('ZPlusJet')
    fakeCategories.append('SameSign')
    fakeCategories.append('MC_truth')

ratioMaker = RatioMaker(inFile, outFile, scale = 19.7)

ratioMaker.get_scale_factors(datasets['QCD2l'] + [bgType], corrected = False)

fakeDict1D = {
    'MuonFakePt':('MuNumerPt', 'MuDenomPt'),
    'MuonFakePtLowJet':('MuNumerPtLowJet', 'MuDenomPtLowJet'),
    'MuonFakePtHighJet':('MuNumerPtHighJet', 'MuDenomPtHighJet'),
    'MuonFakeEta':('MuNumerEta', 'MuDenomEta'),
    'MuonFakeEtaLowJet':('MuNumerEtaLowJet', 'MuDenomEtaLowJet'),
    'MuonFakeEtaHighJet':('MuNumerEtaHighJet', 'MuDenomEtaHighJet'),
    'MuonFakeMet':('MuNumerMet', 'MuDenomMet'),
    'MuonFakeHT':('MuNumerHT', 'MuDenomHT'),
    #'MuonFakeJetMult':('MuNumerJetMult', 'MuDenomJetMult'),
    'ElectronFakePt':('EleNumerPt', 'EleDenomPt'),
    'ElectronFakePtLowJet':('EleNumerPtLowJet', 'EleDenomPtLowJet'),
    'ElectronFakePtHighJet':('EleNumerPtHighJet', 'EleDenomPtHighJet'),
    #'ElectronFakeEta':('EleNumerEta', 'EleDenomEta'),
    'ElectronFakeMet':('EleNumerMet', 'EleDenomMet'),
}

fakeDict2D = {
    'MuonFake':('MuNumer', 'MuDenom'),
    'ElectronFake':('EleNumer', 'EleDenom')
}

for category in fakeCategories:
    print category
    for dataset in datasets[category]:
        ratioMaker.set_category(category)
        ratioMaker.set_dataset(dataset)

        ratioMaker.set_ratio_1D(fakeDict1D)
        ratioMaker.make_1D_ratios(dataset, bgType)

        #ratioMaker.set_ratio_2D(fakeDict2D)
        #ratioMaker.make_2D_ratios(dataset, bgType, doProjections = True)

ratioMaker.write_outfile()

# Combined fakeCategory rates
if doCombination:
    fTest = r.TFile('data/fakeRates_TEST.root', 'UPDATE')
    fTest.mkdir('Combined')
    fTest.cd('Combined')

    #Get 1D histograms
    histList_1D = ['MuonFakePt', 'MuonFakeEta', 'ElectronFakePt', 'ElectronFakeEta']
    histList_2D = ['MuonFake', 'ElectronFake']
    outHists = []
    for hist in histList_1D:
        h1_QCD2l        = fTest.GetDirectory('QCD2l').Get('h1_{0}'.format(hist))
        h1_ZPlusJet     = fTest.GetDirectory('ZPlusJet').Get('h1_{0}'.format(hist))

        h1_QCD2l.SetBit(r.TH1.kIsAverage)    
        h1_ZPlusJet.SetBit(r.TH1.kIsAverage) 

        h1_combined = r.TH1D('h1_{0}'.format(hist), ';p_{T};#varepsilon', h1_QCD2l.GetNbinsX(), h1_QCD2l.GetXaxis().GetXmin(), h1_QCD2l.GetXaxis().GetXmax())
        h1_combined.SetBit(r.TH1.kIsAverage)

        h1_combined.Add(h1_QCD2l)
        h1_combined.Add(h1_ZPlusJet)

        outHists.append(h1_combined)

    for hist in histList_2D:
        h2_QCD2l        = fTest.GetDirectory('QCD2l').Get('h2_{0}'.format(hist))
        h2_ZPlusJet     = fTest.GetDirectory('ZPlusJet').Get('h2_{0}'.format(hist))

        h2_QCD2l.SetBit(r.TH1.kIsAverage)    
        h2_ZPlusJet.SetBit(r.TH1.kIsAverage) 

        h2_combined = r.TH2D('h2_{0}'.format(hist), ';p_{T};#varepsilon', 
                         h2_QCD2l.GetNbinsX(), h2_QCD2l.GetXaxis().GetXmin(), h2_QCD2l.GetXaxis().GetXmax(), 
                         h2_QCD2l.GetNbinsY(), h2_QCD2l.GetYaxis().GetXmin(), h2_QCD2l.GetYaxis().GetXmax())

        h2_combined.SetBit(r.TH1.kIsAverage)

        h2_combined.Add(h2_QCD2l)
        h2_combined.Add(h2_ZPlusJet)

        outHists.append(h2_combined)

        fTest.Write()
        fTest.Close()

