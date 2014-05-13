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
catList = ['inclusive', 'AntiIso3l', 'QCD2l', 'ZPlusJet']
#catList = ['ZPlusJet']

### Samples to be included in stacks ###
samples = {}
samples['inclusive']    = ['PROMPT'] 
samples['ZPlusJet']     = ['PROMPT']
samples['QCD2l']        = ['PROMPT']
samples['AntiIso3l']    = ['PROMPT']
#samples['AntiIso3l']    = ['ZJets', 'ttbarLep', 'ttbarHad', 'WZJets3LNu', 'WbbToLNu'] #'WJetsToLNu']

if doPlots:

    print '\nMaking the plots...\n'

    r.gROOT.SetBatch()

    ### Initialize plot producer ###
    plotter = PlotProducer(inputFile = 'fakeEstimator/histos/{0}.root'.format(batch), savePath = '', scale = LUMIDATA, isAFS = False)
    plotter.set_period(period)
    plotter.set_output_type(plotType)
    plotter.set_save_path('plots/{0}/{1}_{2}_{3}/log'.format(currentDate, selection, batch, suffix))
    plotter.set_clean_fakes(True)

    ### DATASETS ###

    plotter.add_datasets(samples['inclusive'])
    plotter._overlayList.extend(['DATA'])

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
        plotter.add_datasets(samples[category], Clear=True)
        plotter.set_category(category)
        plotter.make_overlays_1D(logScale = doLog, doRatio = doRatio, doEff = doEff)

        plotter.make_overlays_diff([(['DATA_FAKES'], ['MuUnevenPtClosure_{0}'.format(category)]), (['DATA'],['MuNumerPt'])], 'Lepton', 'MuClosurePt') 
        plotter.make_overlays_diff([(['DATA_FAKES'], ['EleUnevenPtClosure_{0}'.format(category)]), (['DATA'],['EleNumerPt'])], 'Lepton', 'EleClosurePt') 

        # Closure plots
        for CR in catList:
            if CR == 'inclusive':
                continue

            plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['MuNumerPt', 'MuUnevenPtClosure_{0}'.format(CR)]), 
                                        (['DATA'],['MuNumerPt'])], 'Lepton', 'MuClosurePt_{0}'.format(CR)) 
            plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['EleNumerPt', 'EleUnevenPtClosure_{0}'.format(CR)]), 
                                        (['DATA'],['EleNumerPt'])], 'Lepton', 'EleClosurePt_{0}'.format(CR)) 
