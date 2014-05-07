#! /usr/bin/env python
import subprocess, shlex, datetime, copy
from multiprocessing import Process

import ROOT as r
from PlotProducer import *
from TableMaker import *

now         = datetime.datetime.now()
currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

### Get command line arguements

if len(sys.argv) > 1:
    batch   = sys.argv[1]
    suffix  = sys.argv[2]
else:
    batch   = '20131018_000148'
    suffix  = 'test'

### This is the config file for manipulating 
### histograms using the PlotProducer class.  

plotType    = '.png'
selection   = 'fakes'

period      = '2012'
LUMIDATA    = 19.712 

doPlots     = True
doLog       = False
doEff       = False
doRatio     = False

### Categories to be plotted ###
catList = ['inclusive', '2l', '3l', 'QCD2l_inclusive', 'ZPlusJet_inclusive']

### Samples to be included in stacks ###
samples = {}
samples['inclusive']    = ['FAKES_2l', 'FAKES_3l'] # ['ZJets', 'ttbar']
samples['ZPlusJet']     = ['FAKES_3l']
samples['QCD2l']        = ['FAKES_2l']
samples['2l']           = ['ZJets', 'ttbarLep', 'ttbarHad', 'Diboson', 'WbbToLNu'] #'WJetsToLNu']
samples['3l']           = ['ZJets', 'ttbarLep', 'ttbarHad', 'WZJets3LNu', 'WbbToLNu'] #'WJetsToLNu']

if doPlots:

    print '\nMaking the plots...\n'

    r.gROOT.SetBatch()

    ### Initialize plot producer ###
    plotter = PlotProducer(inputFile = 'fakeEstimator/histos/{0}.root'.format(batch), savePath = '', scale = LUMIDATA, isAFS = False)
    plotter.set_period(period)
    plotter.set_output_type(plotType)
    plotter.set_save_path('plots/{0}/{1}_{2}_{3}/log'.format(currentDate, selection, batch, suffix))

    ### DATASETS ###

    plotter.add_datasets(samples['ZPlusJet'] + samples['QCD2l'])
    plotter._overlayList.extend(['DATA_FAKES'])
    #plotter._overlayList.extend(['DATA_MUON'])

    plotter.get_scale_factors(addData = [], corrected = False)

    plotter._directoryList1D            = ['Misc', 'Lepton']#, 'CR']
    plotter._variableDict['Misc']       = ['bJetLooseMult', 'bJetMediumMult', 'JetMult', 'Met',
                                           'TagProbeDeltaPhi', 'TagProbePtBalance', 
                                           'MuonIso_1', 'MuonIso_2', 'MuonIso_3',
                                           'MuonIsoRel_Iso', 'MuonIsoRel_AntiIso'
                                           'MuonPt_Iso', 'MuonPt_AntiIso',
                                           'ElectronIso_1', 'ElectronIso_2', 'ElectronIso_3',
                                           'ElectronIsoRel_Iso', 'ElectronIsoRel_AntiIso',
                                           'ElectronPt_Iso', 'ElectronPt_AntiIso'
                                           'TagEleProbeMass', 'TagMuProbeMass']

    plotter._variableDict['Lepton']     = ['MuPassLepPt', 'MuPassLepEta', 
                                           'MuProbeLepPt', 'MuProbeLepEta', 
                                           'MuDenomPt', 'MuDenomEta', 
                                           'MuNumerPt', 'MuNumerEta', 
                                           'MuDenomMet', 'MuNumerMet', 
                                           'MuDenomIsoRel', 'MuNumerIsoRel',
                                           'ElePassLepPt', 'ElePassLepEta', 
                                           'EleProbeLepPt', 'EleProbeLepEta', 
                                           'EleDenomPt', 'EleDenomEta', 
                                           'EleNumerPt', 'EleNumerEta', 
                                           'EleDenomMet', 'EleNumerMet', 
                                           'EleDenomIsoRel', 'EleNumerIsoRel',
                                           'TagLepPt', 'TagLepEta', 'TagIsoRel', 'TagDz', 'TagDxy'] 

    #plotter._variableDict['CR']         = ['MuonIso_1', 'MuonIso_2', 'MuonIso_3', 
    #                                       'MuonIsoRel_Iso', 'MuonIsoRel_AntiIso',
    #                                       'MuonPt_Iso', 'MuonPt_AntiIso',
    #                                       'ElectronIso_1', 'ElectronIso_2', 'ElectronIso_3', 
    #                                       'ElectronIsoRel_Iso', 'ElectronIsoRel_AntiIso',
    #                                       'ElectronPt_Iso', 'ElectronPt_AntiIso']


     ###################   
     ### MAKE PLOTS! ###  
     ###################   

    r.gROOT.SetStyle('Plain')
    r.gStyle.SetOptStat(0)

    ### inclusive ###

    plotter.add_datasets(samples, Clear=True)

    for category in catList:
        plotter.add_datasets(samples[category.split('_', 1)[0]], Clear=True)
        plotter.set_category(category)
        plotter.make_overlays_1D(logScale = doLog, doRatio = doRatio, doEff = doEff)

        # Closure plots
        if category is 'QCD2l_inclusive':
            plotter.make_overlays_diff([(['FAKES_2l', 'DATA_FAKES'], ['MuNumerPt', 'MuUnevenPtClosure']), (['DATA'],['MuNumerPt'])], 'Lepton', 'MuClosurePt') 
            plotter.make_overlays_diff([(['FAKES_2l', 'DATA_FAKES'], ['EleNumerPt', 'EleUnevenPtClosure']), (['DATA'],['EleNumerPt'])], 'Lepton', 'EleClosurePt') 
        if category is 'ZPlusJet_inclusive':
            plotter.make_overlays_diff([(['PROMPT_3l', 'DATA_FAKES'], ['MuNumerPt', 'MuUnevenPtClosure']), (['PASS'],['MuNumerPt'])], 'Lepton', 'MuClosurePt') 
            plotter.make_overlays_diff([(['PROMPT_3l', 'DATA_FAKES'], ['EleNumerPt', 'EleUnevenPtClosure']), (['PASS'],['EleNumerPt'])], 'Lepton', 'EleClosurePt') 

        # X-check plots
        if category in ['2l', '3l']:
            plotter.make_overlays_diff([(['PROMPT_2l', 'FAKEABLE'], ['MuonPt_QCD2l_tight', 'MuonPt_QCD2l_weight']), (['DATA'],['MuonPt_QCD2l_tight'])], 'Lepton', 'MuonFakePt') 
            plotter.make_overlays_diff([(['PROMPT_2l', 'FAKEABLE'], ['MuonEta_QCD2l_tight', 'MuonEta_QCD2l_weight']), (['DATA'],['MuonEta_QCD2l_tight'])], 'Lepton', 'MuonFakeEta') 
            plotter.make_overlays_diff([(['PROMPT_2l', 'FAKEABLE'], ['MuonMet_QCD2l_tight', 'MuonMet_QCD2l_weight']), (['DATA'],['MuonMet_QCD2l_tight'])], 'Lepton', 'MuonFakeMet') 
            plotter.make_overlays_diff([(['PROMPT_2l', 'FAKEABLE'], ['MuonJetMult_QCD2l_tight', 'MuonJetMult_QCD2l_weight']), (['DATA'],['MuonJetMult_QCD2l_tight'])], 'Lepton', 'MuonFakeJetMult') 

            plotter.make_overlays_diff([(['PROMPT_2l', 'FAKEABLE'], ['ElectronPt_QCD2l_tight', 'ElectronPt_QCD2l_weight']), (['DATA'],['ElectronPt_QCD2l_tight'])], 'Lepton', 'ElectronFakePt') 
            plotter.make_overlays_diff([(['PROMPT_2l', 'FAKEABLE'], ['ElectronEta_QCD2l_tight', 'ElectronEta_QCD2l_weight']), (['DATA'],['ElectronEta_QCD2l_tight'])], 'Lepton', 'ElectronFakeEta') 
            plotter.make_overlays_diff([(['PROMPT_2l', 'FAKEABLE'], ['ElectronMet_QCD2l_tight', 'ElectronMet_QCD2l_weight']), (['DATA'],['ElectronMet_QCD2l_tight'])], 'Lepton', 'ElectronFakeMet') 
            plotter.make_overlays_diff([(['PROMPT_2l', 'FAKEABLE'], ['ElectronJetMult_QCD2l_tight', 'ElectronJetMult_QCD2l_weight']), (['DATA'],['ElectronJetMult_QCD2l_tight'])], 'Lepton', 'ElectronFakeJetMul') 

