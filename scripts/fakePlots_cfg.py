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
    batch = sys.argv[1]
else:
    batch = '20131018_000148'

### This is the config file for manipulating 
### histograms using the PlotProducer class.  

plotType    = '.png'
selection   = 'fakes'

period      = '2012'
LUMIDATA    = 19.712 
doLog       = True

doPlots     = True
doRates     = False

### Categories to be plotted ###
catList = ['inclusive', 'low_met', 'high_met']
#catList = ['low_met']

### Samples to be included in stacks ###
samples = ['ZJets', 'ttbar']

if doPlots:

    print '\nMaking the plots...\n'

    r.gROOT.SetBatch()

    ### Initialize plot producer ###
    plotter = PlotProducer(inputFile = 'fakeEstimator/histos/{0}.root'.format(batch), savePath = '', scale = LUMIDATA, isAFS = False)
    plotter.set_period(period)
    plotter.set_output_type(plotType)
    plotter.set_save_path('plots/{0}/{1}_{2}/log'.format(currentDate, selection, batch))

    ### DATASETS ###

    plotter.add_datasets(samples)
    plotter._overlayList.extend(['DATA'])

    plotter.get_scale_factors(corrected = False)

    plotter._directoryList1D            = ['Misc', 'Lepton']
    plotter._variableDict['Misc']       = ['bJetMult', 'JetMult', 'TagProbeDeltaPhi', 'TagProbePtBalance', 'Met']
    plotter._variableDict['Lepton']     = ['MuPassLepPt', 'MuPassLepEta', 
                                           'MuProbeLepPt', 'MuProbeLepEta', 
                                           'MuDenomPt', 'MuDenomEta', 
                                           'MuNumerPt', 'MuNumerEta', 
                                           'TagLepPt', 'TagLepEta' ] 


     ###################   
     ### MAKE PLOTS! ###  
     ###################   

    r.gROOT.SetStyle('Plain')
    r.gStyle.SetOptStat(0)

    ### inclusive ###

    plotter.add_datasets(samples, Clear=True)

    for category in catList:
        plotter.set_category(category)
        plotter.make_overlays_1D(logScale = doLog, doRatio = True, doEff = False)

