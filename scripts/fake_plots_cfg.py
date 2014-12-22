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

### Samples to be included in stacks ###
samples = {}
samples['inclusive']    = ['Rare', 'ZJets', 'ttbar', 'WJets', 'WZJets3LNu', 'QCD']
samples['ZPlusJet']     = ['Rare', 'WZJets3LNu', 'ZJets', 'WJets']
samples['QCD2l']        = ['ZJets', 'WJets', 'QCD', 'ttbar']
samples['SameSign']     = ['Rare', 'WZJets3LNu', 'ttbar', 'ZJets', 'WJets', 'QCD']
#samples['SameSign']     = ['WZJets3LNu', 'ttbar', 'ZJets', 'W1JetsToLNu', 'W2JetsToLNu', 'W3JetsToLNu', 'W4JetsToLNu']

if doPlots:

    print '\nMaking the plots...\n'

    r.gROOT.SetBatch()

    ### Initialize plot producer ###
    plotter = PlotProducer(inputFile = 'fakeEstimator/histos/{0}.root'.format(batch), savePath = '', scale = LUMIDATA, isAFS = False)
    plotter.set_period(period)
    plotter.set_output_type(plotType)
    plotter.set_save_path('plots/{0}/{1}_{2}_{3}/log'.format(currentDate, selection, batch, suffix))

    ### DATASETS ###

    plotter.add_datasets(samples['inclusive'])
    plotter._overlayList.extend(['DATA'])
    plotter.get_scale_factors(corrected = False)

    plotter._directoryList1D            = ['Muon', 'Electron']

    plotter._variableDict['Muon']       = ['MuPassLepPt', 'MuPassLepEta', 'MuPassTransverseMass', 
                                           'MuProbeLepPt', 'MuProbeLepEta', 'MuProbeTransverseMass',
                                           'MuNumerPt', 'MuNumerEta', 'MuNumerMet', 'MuNumerIsoRel',
                                           'MuDenomPt', 'MuDenomEta', 'MuDenomMet', 'MuDenomIsoRel',
                                           'MuNumerPtLowJet', 'MuNumerPtHighJet',
                                           'MuDenomPtLowJet', 'MuDenomPtHighJet',
                                           'MuDenomIsoRelBin1', 'MuDenomIsoRelBin2',
                                           'MuFailLepPt', 'MuFailLepEta', 'MuFailTransverseMass', 
                                           'MuTagLepPt', 'MuTagLepEta', 'MuTagIsoRel', 'MuTagDz', 'MuTagDxy', 
                                           'TagMuPassMass', 'TagMuPassDeltaR',
                                           'TagMuProbeMass', 'TagMuProbeDeltaR',
                                           'TagMuProbeDeltaPhi', 'TagMuProbePtBalance', 
                                           'TagMuFailMass', 'TagMuFailDeltaR',
                                           'MuJetMult_tight', 'MuJetMult_inclusive', 'MuJetMult_fail']

    plotter._variableDict['Electron']   = ['ElePassLepPt', 'ElePassLepEta', 'ElePassTransverseMass', 
                                           'EleProbeLepPt', 'EleProbeLepEta', 'EleProbeTransverseMass',
                                           'EleDenomPt', 'EleDenomEta', 
                                           'EleNumerPt', 'EleNumerEta', 
                                           'EleDenomMet', 'EleNumerMet', 
                                           'EleDenomIsoRel', 'EleNumerIsoRel',
                                           'EleDenomIsoRelBin1', 'EleDenomIsoRelBin2',
                                           'TagEleProbeMass',
                                           'TagEleProbeDeltaPhi', 'TagEleProbePtBalance', 
                                           'EleTagLepPt', 'EleTagLepEta', 'EleTagIsoRel', 'EleTagDz', 'EleTagDxy'
                                           ]

    plotter._variableDict['GenStudies'] = ['MatchedMuJetFlavor_inclusive', 'MatchedMuJetBDiscr_inclusive',
                                           'LightMatchedMuIsoRel_inclusive', 'LightMatchedMuPt_inclusive', 
                                           'LightMatchedMuEta_inclusive', 'LightMatchedMuMet_inclusive', 'LightMatchedMuMT_inclusive',
                                           'HeavyMatchedMuIsoRel_inclusive', 'HeavyMatchedMuPt_inclusive', 'HeavyMatchedMuEta_inclusive', 
                                           'HeavyMatchedMuMet_inclusive', 'HeavyMatchedMuMT_inclusive',
                                           'MatchedMuJetFlavor_tight', 'MatchedMuJetBDiscr_tight',
                                           'LightMatchedMuIsoRel_tight', 'LightMatchedMuPt_tight', 
                                           'LightMatchedMuEta_tight', 'LightMatchedMuMet_tight', 'LightMatchedMuMT_tight',
                                           'HeavyMatchedMuIsoRel_tight', 'HeavyMatchedMuPt_tight', 'HeavyMatchedMuEta_tight', 
                                           'HeavyMatchedMuMet_tight', 'HeavyMatchedMuMT_tight',
                                           'MatchedMuJetFlavor_fail', 'MatchedMuJetBDiscr_fail',
                                           'LightMatchedMuIsoRel_fail', 'LightMatchedMuPt_fail', 
                                           'LightMatchedMuEta_fail', 'LightMatchedMuMet_fail', 'LightMatchedMuMT_fail',
                                           'HeavyMatchedMuIsoRel_fail', 'HeavyMatchedMuPt_fail', 'HeavyMatchedMuEta_fail', 
                                           'HeavyMatchedMuMet_fail', 'HeavyMatchedMuMT_fail']
                                           


     ###################   
     ### MAKE PLOTS! ###  
     ###################   

    r.gROOT.ProcessLine('.L ./scripts/tdrStyle.C')
    r.setTDRStyle()

    #r.TGaxis.SetMaxDigits(3)
    r.gStyle.SetOptStat(0)

    ### Categories to be plotted ###
    catList = ['QCD2l', 'ZPlusJet', 'SameSign']

    for category in catList:
        plotter._directoryList1D = ['Muon', 'Electron', 'GenStudies']
        plotter.add_datasets(samples[category], Clear=True)
        plotter.set_category(category)
        plotter.make_overlays_1D(logScale = doLog, doRatio = doRatio, doEff = doEff)

        # Closure plots
        #for CR in ['QCD2l', 'ZPlusJet']:

        #    plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['MuNumerEta', 'MuUnevenEtaClosure_{0}'.format(CR)]), 
        #                                (['DATA'],['MuNumerEta'])], 'Muon', 'MuClosureEta_{0}'.format(CR)) 
        #    plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['EleNumerEta', 'EleUnevenEtaClosure_{0}'.format(CR)]), 
        #                                (['DATA'],['EleNumerEta'])], 'Electron', 'EleClosureEta_{0}'.format(CR)) 

        #    plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['MuNumerPt', 'MuUnevenPtClosure_{0}'.format(CR)]), 
        #                                (['DATA'],['MuNumerPt'])], 'Muon', 'MuClosurePt_{0}'.format(CR)) 
        #    plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['EleNumerPt', 'EleUnevenPtClosure_{0}'.format(CR)]), 
        #                                (['DATA'],['EleNumerPt'])], 'Electron', 'EleClosurePt_{0}'.format(CR)) 

        #    plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['MuNumerMet', 'MuMetClosure_{0}'.format(CR)]), 
        #                                (['DATA'],['MuNumerMet'])], 'Muon', 'MuClosureMet_{0}'.format(CR)) 
        #    plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['EleNumerMet', 'EleMetClosure_{0}'.format(CR)]), 
        #                                (['DATA'],['EleNumerMet'])], 'Electron', 'EleClosureMet_{0}'.format(CR)) 

        #    plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['MuNumerJetMult', 'MuJetMultClosure_{0}'.format(CR)]), 
        #                                (['DATA'],['MuNumerJetMult'])], 'Muon', 'MuClosureJetMult_{0}'.format(CR)) 
        #    plotter.make_overlays_diff([(['PROMPT', 'DATA_FAKES'], ['EleNumerJetMult', 'EleJetMultClosure_{0}'.format(CR)]), 
        #                                (['DATA'],['EleNumerJetMult'])], 'Electron', 'EleClosureJetMult_{0}'.format(CR)) 

        #plotter._directoryList1D = ['GenStudies']
        #plotter.make_stacks_by_category(logScale = doLog)
