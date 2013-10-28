#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user'
executable  = 'execBatch.csh'

selection   = 'fcnc'
period      = '2012'
doData      = False
doBG        = False
doSignal    = False
doFakes     = False

# Config from command line #

if len(sys.argv) > 0:
    period = sys.argv[1]
if 'data' in sys.argv[1:]:
    doData = True
if 'bg' in sys.argv[1:]:
    doBG = True
if 'signal' in sys.argv[1:]:
    doSignal = True
if 'fakes' in sys.argv[1:]:
    doFakes = True

############################

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

data    = []
fakes   = []
bg      = []
signal  = []


if period == '2012':
    data.extend([
        cfg('muon_2012A', dCache+'/naodell/nuTuples_v7_4/DoubleMu_Run2012A', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', dCache+'/naodell/nuTuples_v7_4/DoubleMu_Run2012B', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', dCache+'/naodell/nuTuples_v7_4/DoubleMu_Run2012C', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', dCache+'/naodell/nuTuples_v7_4/DoubleMu_Run2012D', 15, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', dCache+'/naodell/nuTuples_v7_4/DoubleElectron_Run2012A', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', dCache+'/naodell/nuTuples_v7_4/DoubleElectron_Run2012B', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', dCache+'/naodell/nuTuples_v7_4/DoubleElectron_Run2012C', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', dCache+'/naodell/nuTuples_v7_4/DoubleElectron_Run2012D', 15, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', dCache+'/naodell/nuTuples_v7_4/MuEG_Run2012A', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', dCache+'/naodell/nuTuples_v7_4/MuEG_Run2012B', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', dCache+'/naodell/nuTuples_v7_4/MuEG_Run2012C', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', dCache+'/naodell/nuTuples_v7_4/MuEG_Run2012D', 15, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        #cfg('ZJets_M-50', dCache+'/naodell/nuTuples_v7_4/DYJets', 30, 'ZJets_M-50 muon 2012'),
        #cfg('ZJets_M-10To50', dCache+'/naodell/nuTuples_v7_4/DYJets_M-10To50', 10, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZbbToLL', dCache+'/naodell/nuTuples_v7_4/ZbbToLL', 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', dCache+'/naodell/nuTuples_v7_4/ZGToLLG', 10, 'ZG muon 2012'),
        #cfg('WJets', dCache+'/naodell/nuTuples_v7_4/WJetsToLNu', 20, 'WJets muon 2012'),
        #cfg('WbbToLNu', dCache+'/naodell/nuTuples_v7_4/WbbJetsToLNu', 20, 'WbbToLNu muon 2012'),
        #cfg('WG', dCache+'/naodell/nuTuples_v7_4/WGToLNuG', 10, 'WG muon 2012'),

        #cfg('ttbar', dCache+'/naodell/nuTuples_v7_4/TTJets', 30, 'ttbar muon 2012'),
        #cfg('tbarW', dCache+'/naodell/nuTuples_v7_4/Tbar_tW', 5, 'tbarW muon 2012'),
        #cfg('tW', dCache+'/naodell/nuTuples_v7_4/T_tW', 5, 'tW muon 2012'),
        #cfg('WWJets2L2Nu', dCache+'/naodell/nuTuples_v7_4/WWJetsTo2L2Nu', 5, 'WWJets2L2Nu muon 2012'),
        #cfg('ZZJets2L2Nu', dCache+'/naodell/nuTuples_v7_4/ZZJetsTo2L2Nu', 5, 'ZZJets2L2Nu muon 2012'),
        cfg('ZZJets2L2Q', dCache+'/naodell/nuTuples_v7_4/ZZJetsTo2L2Q', 5, 'ZZJets2L2Q muon 2012')

        #cfg('QCD_20-30_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_20_30_EMEnriched', 10, 'QCD_20-30_EM muon 2012'),
        #cfg('QCD_30-80_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_30_80_EMEnriched', 10, 'QCD_30-80_EM muon 2012'),
        #cfg('QCD_80-170_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_80_170_EMEnriched', 10, 'QCD_80-170_EM muon 2012'),
        #cfg('QCD_170-250_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_170_250_EMEnriched', 10, 'QCD_170-250_EM muon 2012'),
        #cfg('QCD_250-350_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_250_350_EMEnriched', 10, 'QCD_250-350_EM muon 2012'),
        #cfg('QCD_350_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_350_EMEnriched', 10, 'QCD_350_EM muon 2012'),
        #cfg('QCD_20_MU', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_20_MuEnrichedPt_15', 10, 'QCD_20_MU muon 2012')

        ])

    signal.extend([
        cfg('FCNC_M125_tHj', dCache+'/naodell/nuTuples_v6_8TeV/FCNH_M125_t', 1, 'FCNC_M125_t mc 2012'),
        cfg('FCNC_M125_tbarHj', dCache+'/naodell/nuTuples_v6_8TeV/FCNH_M125_tbar', 1, 'FCNC_M125_tbar mc 2012')
        ])


inputSamples = []

if doData:
    inputSamples.extend(data)
if doBG:
    inputSamples.extend(bg)
if doSignal:
    inputSamples.extend(signal)
if doFakes:
    inputSamples.extend(fakes)

if len(inputSamples) is not 0:
    batcher = b.BatchMaster(inputSamples, shortQueue = False, stageDir = 'batchStage', executable = executable, selection = selection + '_' + period)
    batcher.submit_to_batch()

