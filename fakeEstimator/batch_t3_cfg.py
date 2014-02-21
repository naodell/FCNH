#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
mcDir       = '/tthome/naodell/storage/data/nuTuples_v9.6_8TeV/MC'
dataDir     = '/tthome/bpollack/storage/nuTuples_v9.6_8TeV/Data/'
executable  = 'execBatch.sh'

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

############################

''' 
    Set job configurations.  The order of arguments is:
    (Dataset, path to data, number of jobs, arguments to pass to executable, output directory name)
'''

data    = []
bg      = []
signal  = []


if period == '2012':
    data.extend([
        cfg('muon_2012A', '{0}/DoubleMu_Run2012A'.format(dataDir), 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', '{0}/DoubleMu_Run2012B'.format(dataDir), 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', '{0}/DoubleMu_Run2012C'.format(dataDir), 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', '{0}/DoubleMu_Run2012D'.format(dataDir), 15, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', '{0}/DoubleElectron_Run2012A'.format(dataDir), 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', '{0}/DoubleElectron_Run2012B'.format(dataDir), 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', '{0}/DoubleElectron_Run2012C'.format(dataDir), 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', '{0}/DoubleElectron_Run2012D'.format(dataDir), 15, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', '{0}/MuEG_Run2012A'.format(dataDir), 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', '{0}/MuEG_Run2012B'.format(dataDir), 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', '{0}/MuEG_Run2012C'.format(dataDir), 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', '{0}/MuEG_Run2012D'.format(dataDir), 15, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        cfg('ZJets_M-50', mcDir+'/DYJetsToLL_M-50', 50, 'ZJets_M-50 muon 2012'),
        cfg('ZJets_M-10To50', mcDir+'/DYJetsToLL_M-10To50filter', 35, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZbbToLL', mcDir+'/ZbbToLL', 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', mcDir+'/ZGToLLG', 10, 'ZG muon 2012'),
        #cfg('WJets', mcDir+'/WJetsToLNu', 20, 'WJets muon 2012'),
        #cfg('WGStarLNu2E', mcDir+'/WGStarToLNu2E', 5, 'WGStarLNu2E muon 2012'),
        #cfg('WGStarLNu2Mu', mcDir+'/WGStarToLNu2Mu', 5, 'WGStarLNu2Mu muon 2012'),
        #cfg('WGStarLNu2Tau', mcDir+'/WGStarToLNu2Tau', 5, 'WGStarLNu2Tau muon 2012'),
        #cfg('WbbToLNu', mcDir+'/WbbJetsToLNu', 20, 'WbbToLNu muon 2012'),
        #cfg('WbbToLNu', mcDir+'/WbbToLL', 20, 'WbbToLNu muon 2012'),
        #cfg('WG', mcDir+'/WGToLNuG', 10, 'WG muon 2012'),

        cfg('ttbarHad', mcDir+'/TTJets', 30, 'ttbarHad muon 2012'),
        cfg('ttbarLep', mcDir+'/TTJets', 30, 'ttbarLep muon 2012'),
        cfg('tbarW', mcDir+'/Tbar_tW', 5, 'tbarW muon 2012'),
        cfg('tW', mcDir+'/T_tW', 5, 'tW muon 2012'),
        cfg('t_t-channel', mcDir+'/T_t', 5, 't_t-channel muon 2012'),
        cfg('tbar_t-channel', mcDir+'/Tbar_t', 5, 'tbar_t-channel muon 2012'),

        cfg('ZZJets2L2Nu', mcDir+'/ZZJetsTo2L2Nu', 5, 'ZZJets2L2Nu muon 2012'),
        #cfg('ZZJets2L2Q', mcDir+'/ZZJetsTo2L2Q', 5, 'ZZJets2L2Q muon 2012'),
        #cfg('ZZJets4L', mcDir+'/ZZJetsTo4L', 5, 'ZZJets4L muon 2012'),
        cfg('ZZTo4e', mcDir+'/ZZTo4e', 5, 'ZZ4e muon 2012'),
        cfg('ZZTo4mu', mcDir+'/ZZTo4mu', 5, 'ZZ4mu muon 2012'),
        cfg('ZZTo4tau', mcDir+'/ZZTo4tau', 5, 'ZZ4tau muon 2012'),
        cfg('ZZTo2e2mu', mcDir+'/ZZTo2e2mu', 5, 'ZZ2e2mu muon 2012'),
        cfg('ZZTo2e2tau', mcDir+'/ZZTo2e2tau', 5, 'ZZ2e2tau muon 2012'),
        cfg('ZZTo2mu2tau', mcDir+'/ZZTo2mu2tau', 5, 'ZZ2mu2tau muon 2012'),
        cfg('WWJets2L2Nu', mcDir+'/WWJetsTo2L2Nu', 5, 'WWJets2L2Nu muon 2012'),
        cfg('WZJets3LNu', mcDir+'/WZJetsTo3LNu', 5, 'WZJets3LNu muon 2012'),
        cfg('WZJets2L2Q', mcDir+'/WZJetsTo2L2Q', 10, 'WZJets2L2Q muon 2012'),
        ])


inputSamples = []

if doData:
    inputSamples.extend(data)
if doBG:
    inputSamples.extend(bg)

if len(inputSamples) is not 0:
    batcher = b.BatchMaster(inputSamples, shortQueue = False, stageDir = 'batchStage', executable = executable, selection = selection + '_' + period)
    batcher.submit_to_batch()

