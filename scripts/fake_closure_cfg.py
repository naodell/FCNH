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
    print 'Need to specify input batch and suffix!!!  Exiting...'
    exit()

### This is the config file for manipulating 
### histograms using the PlotProducer class.  

plotType    = '.png'
selection   = 'fcnh'

cutList     = ['1_preselection']
#cutList.extend(['2_Z_veto', '3_jet', '4_MET'])
fakeType    = 'eFakes'

period      = '2012'
LUMIDATA    = 19.712 

doPlots     = True
doLog       = True
doEff       = False
doRatio     = False
do1D        = True
do2D        = False

doYields    = False

doSS        = True
do3l        = True

### Categories to be plotted ###
catSS       = ['ss_inclusive']
catSS.extend(['ss_mumu', 'ss_ee', 'ss_emu'])
cat3l       = ['3l_inclusive']
cat3l.extend(['3l_OSSF', '3l_SSSF'])
cat3l.extend(['3l_eee', '3l_eemu', '3l_emumu', '3l_mumumu'])

### Samples to be included in stacks ###
samples     = {'all':[], 'inclusive':[], '3l':[], 'ss':[]}

#samples['all'].append('higgs')
samples['all'].append('Triboson')
samples['all'].append('ttV')
samples['all'].append('Diboson')
samples['all'].append('top')
samples['all'].append('ZJets')
#samples['all'].append('WJets')
#samples['all'].append('QCD')
#samples['all'].extend(['ZbbToLL', 'WbbToLNu']) #, 'ZGstar'])

#samples['3l'].append(fakeType + '_Triboson')
samples['3l'].append(fakeType + '_ttV')
samples['3l'].append(fakeType + '_ZZ4l')
samples['3l'].append(fakeType + '_WZJets3LNu')
samples['3l'].append(fakeType + '_ttbarHad')
samples['3l'].append(fakeType + '_ttbarLep')
samples['3l'].append(fakeType + '_ZJets')

#samples['ss'].append(fakeType + '_WWJets2L2Nu')
#samples['ss'].append(fakeType + '_ZZJets2L2Nu')
#samples['ss'].append(fakeType + '_ZZJets2L2Q')
samples['ss'].append(fakeType + '_Diboson')
samples['ss'].append(fakeType + '_ttbarLep')
samples['ss'].append(fakeType + '_ttbarHad')
samples['ss'].append(fakeType + '_ZJets')
#samples['ss'].append(fakeType + '_WJets')
#samples['ss'].append(fakeType + '_QCD')

p_plot = []

if doPlots:

    print '\nMaking the plots...\n'

    r.gROOT.SetBatch()

    ### Initialize plot producer ###
    plotter = PlotProducer(inputFile = 'fcncAnalysis/combined_histos/{0}_cut1_{1}_{2}.root'.format(selection, period, batch), savePath = '', scale = LUMIDATA, isAFS = False)
    plotter.set_period(period)
    plotter.set_output_type(plotType)
    plotter.set_clean_fakes(False)

    ### DATASETS ###
    ### Specify the datasets you wish to stack 
    ### and overlay accordingly. 

    plotter.add_datasets(samples['all'])
    plotter.get_scale_factors()

    ### VARIABLES ###
    ### First specify the directories in which your
    ### histograms are stored.  If directories are 
    ### not used enter '' as the only entry.  Then 
    ### list all of the variable names you wish to 
    ### plot while giving a key value which is the 
    ### directory that they are located in as a key.

    plotter._directoryList1D            = ['Misc', 'Lepton', 'Lep+Jet', 'Dilepton', 'DileptonOS', 'Trilepton', 'MET', 'Jet', 'GEN', '4l']

    plotter._variableDict['Misc']       = ['PvMult', 'YieldByCut', 'YieldByCutRaw', 'EventWeight', 'TriggerStatus', 
                                            'FakeWeightUncertainty', 'BDT', 'FakeCategory']

    plotter._variableDict['Lepton']     = ['LeptonCharge', 'LeptonFlavor', 
                                           'Lepton1Pt', 'Lepton2Pt','Lepton3Pt',
                                           'Lepton1Eta', 'Lepton2Eta', 'Lepton3Eta',
                                           'ElectronPt', 'ElectronEta',
                                           'MuonPt', 'MuonEta',
                                           'Lepton1dxy', 'Lepton1dz',
                                           'Lepton2dxy', 'Lepton2dz',
                                           'Lepton3dxy', 'Lepton3dz',
                                           'LeptonMult', 'fakeableOverlapMult'
                                           'OverlapEleMu', 'MuEleDeltaR']
                                           #'Lepton1Phi', 'Lepton2Phi', 'Lepton3Phi']

    plotter._variableDict['Dilepton']   = ['DileptonMass21', 'DileptonTransMass21', 'DileptonQt21',
                                           'DileptonDeltaPhi21', 'DileptonDeltaEta21', 'DileptonDeltaR21', 'DileptonDeltaPt21',
                                           'DileptonMass31', 'DileptonTransMass31', 'DileptonQt31', 
                                           'DileptonDeltaPhi31', 'DileptonDeltaEta31', 'DileptonDeltaR31', 'DileptonDeltaPt31',
                                           'DileptonMass32', 'DileptonTransMass32', 'DileptonQt32', 
                                           'DileptonDeltaPhi32', 'DileptonDeltaEta32', 'DileptonDeltaR32', 'DileptonDeltaPt32']

    plotter._variableDict['DileptonOS'] = ['DileptonOSMass', 'DileptonOSTransMass', 'DileptonOSBalance',
                                           'DileptonOSQt', 'DileptonOSDeltaPt', 'DileptonOSDeltaR', 
                                           'DileptonOSDeltaEta', 'DileptonOSDeltaPhi'] 

    plotter._variableDict['Trilepton']  = ['DileptonLepDeltaR', 'DileptonLepDeltaPhi', 'DileptonLepDeltaEta', 
                                           'Lep3MetMT', 'TrileptonMass', 'TrileptonPt']

    plotter._variableDict['Lep+Jet']    = ['Lepton1BJetDeltaPhi', 'Lepton1BJetDeltaEta', 'Lepton1BJetDeltaR', 'Lepton1BJetDeltaPt',
                                           'Lepton2BJetDeltaPhi', 'Lepton2BJetDeltaEta', 'Lepton2BJetDeltaR', 'Lepton2BJetDeltaPt',
                                           'Lepton3BJetDeltaPhi', 'Lepton3BJetDeltaEta', 'Lepton3BJetDeltaR', 'Lepton3BJetDeltaPt',
                                           'Lepton1JetDeltaPhi', 'Lepton1JetDeltaEta', 'Lepton1JetDeltaR', 'Lepton1JetDeltaPt',
                                           'Lepton2JetDeltaPhi', 'Lepton2JetDeltaEta', 'Lepton2JetDeltaR', 'Lepton2JetDeltaPt',
                                           'Lepton3JetDeltaPhi', 'Lepton3JetDeltaEta', 'Lepton3JetDeltaR', 'Lepton3JetDeltaPt',
                                           'DileptonJetDeltaPhi', 'DileptonJetDeltaEta', 'DileptonJetDeltaR', 'DileptonJetDeltaPt',
                                           'DileptonBJetDeltaPhi', 'DileptonBJetDeltaEta', 'DileptonBJetDeltaR', 'DileptonBJetDeltaPt',
                                           'OverlapJetMult'
                                          ]

    plotter._variableDict['Jet']        = ['Jet1Pt', 'Jet2Pt',# 'Jet3Pt',
                                           'Jet1Eta', 'Jet2Eta',# 'Jet3Eta',
                                           #'Jet1Phi', 'Jet2Phi', 'Jet3Phi',
                                           'BJet1BDiscr', 'BJet1Pt', 'BJet1Eta', #'BJet1Phi', 
                                           'BJet2BDiscr', 'BJet2Pt', 'BJet2Eta', #'BJet2Phi',
                                           'HT', 'HTs', 'EventBalance', 'Centrality',
                                           'JetBJetDeltaPhi', 'JetBJetDeltaEta', 'JetBJetDeltaR',
                                           'JetMultCharge', 'JetMult', 'BJetMult', 'AllJetMult',
                                           'MatchedMuJetBDiscr', 'MatchedEleJetBDiscr']

    plotter._variableDict['MET']        = ['Met', 'MHT', 'METLD', 'MHT-MET', 'MetPhi', 'MetSumEt',
                                           'MetLepton1DeltaPhi', 'MetLepton2DeltaPhi', 'MetLepton3DeltaPhi'
                                           'MetLepDeltaPhiMin', 'nearLepIndex', 'ProjectedMet',
                                           'MetFakeableDeltaPhi'] 

    plotter._variableDict['GEN']        = ['GenChargeMisId', 'GenMisIdPt', 'GenMisIdEta',
                                           'GenDeltaR', 'GenBalance']

    plotter._variableDict['4l']         = ['4lMass', '4lPt', '4lSumPt', '4lMet']



     ###################   
     ### MAKE PLOTS! ###  
     ###################   

    r.gROOT.SetStyle('Plain')
    r.gStyle.SetOptStat(0)

    ### 3l selection ###
    if do3l:

        plotter_3l = copy.deepcopy(plotter)
        plotter_3l.add_datasets(samples['3l'], Clear=True)
        plotter_3l._overlayList = [fakeType]

        for i, cut in enumerate(cutList):
            inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

            if doLog:
                outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
            else:
                outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

            plotter_3l.make_save_path(outFile, clean=True)

            for category in cat3l:
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(plotter_3l, category, inFile, outFile, do1D, do2D, False, doLog, doRatio, doEff)))

    ### ss selection ###
    if doSS:
        ss_plotter = copy.deepcopy(plotter)
        ss_plotter.add_datasets(samples['ss'], Clear=True)
        ss_plotter._overlayList = [fakeType]

        for i, cut in enumerate(cutList):
            inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

            if doLog:
                outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
            else:
                outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

            ss_plotter.make_save_path(outFile, clean=True)

            for category in catSS:
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(ss_plotter, category, inFile, outFile, do1D, do2D, False, doLog, doRatio, doEff)))

    ### ZPlusFake control region
    ZFake_plotter = copy.deepcopy(plotter)
    ZFake_plotter.add_datasets(samples['3l'],  Clear=True)
    ZFake_plotter._overlayList = [fakeType]

    inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 10, period, batch)

    if doLog:
        outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'CR_ZFake')
    else:
        outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'CR_ZFake')

    ZFake_plotter.make_save_path(outFile, clean=True)

    for category in cat3l:
        p_plot.append(Process(name = 'CR_ZFake/' + category, target = plotter_wrapper, args=(ZFake_plotter, category, inFile, outFile, do1D, False, False, doLog, doRatio, False)))

### End of configuration for PlotProducer ###

for process in p_plot:
    print 'Plotting {0}'.format(process.name)
    process.start()

for process in p_plot:
    process.join()

print '\n'

