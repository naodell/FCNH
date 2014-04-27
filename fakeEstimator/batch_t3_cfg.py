#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
dataDir     = '/tthome/bpollack/storage/nuTuples_v9.6_8TeV/Data'
mcDir       = '/tthome/naodell/storage/data/nuTuples_v9.6_8TeV/MC'
executable  = 'execBatch.sh'

selection   = 'fcnc'
period      = '2012'
doData      = False
doBG        = False
doSignal    = False
doFakes     = False
mcTrigger   = 'muon'

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
        cfg('muon_2012A', '{0}/DoubleMu_Run2012A'.format(dataDir), 20, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', '{0}/DoubleMu_Run2012B'.format(dataDir), 20, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', '{0}/DoubleMu_Run2012C'.format(dataDir), 20, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', '{0}/DoubleMu_Run2012D'.format(dataDir), 25, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', '{0}/DoubleElectron_Run2012A'.format(dataDir), 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', '{0}/DoubleElectron_Run2012B'.format(dataDir), 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', '{0}/DoubleElectron_Run2012C'.format(dataDir), 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', '{0}/DoubleElectron_Run2012D'.format(dataDir), 25, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', '{0}/MuEG_Run2012A'.format(dataDir), 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', '{0}/MuEG_Run2012B'.format(dataDir), 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', '{0}/MuEG_Run2012C'.format(dataDir), 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', '{0}/MuEG_Run2012D'.format(dataDir), 25, 'DATA_MUEG muEG 2012')
    ])

    bg.extend([
        cfg('ZJets_M-50',      '{0}/DYJetsToLL_M-50'.format(mcDir),            30,  'ZJets_M-50      {0}  2012'.format(mcTrigger)),
        cfg('ZJets_M-10To50',  '{0}/DYJetsToLL_M-10To50filter'.format(mcDir),  20,  'ZJets_M-10To50  {0}  2012'.format(mcTrigger)),
        cfg('WbbToLNu',        '{0}/WbbJetsToLNu'.format(mcDir),     20,  'WbbToLNu        {0}  2012'.format(mcTrigger)),
        cfg('WjetToLNu',       '{0}/WJetsToLNu'.format(mcDir),       20,  'WJetsToLNu      {0}  2012'.format(mcTrigger)),

        cfg('ttbarHad',        '{0}/TTJets'.format(mcDir),           30,  'ttbarHad        {0}  2012'.format(mcTrigger)),
        cfg('ttbarLep',        '{0}/TTJets'.format(mcDir),           30,  'ttbarLep        {0}  2012'.format(mcTrigger)),
        cfg('tbarW',           '{0}/Tbar_tW'.format(mcDir),          5,   'tbarW           {0}  2012'.format(mcTrigger)),
        cfg('tW',              '{0}/T_tW'.format(mcDir),             5,   'tW              {0}  2012'.format(mcTrigger)),

        cfg('WZJets3LNu',      '{0}/WZJetsTo3LNu'.format(mcDir),     5,   'WZJets3LNu      {0}  2012'.format(mcTrigger)),
        cfg('WWJets2L2Nu',     '{0}/WWJetsTo2L2Nu'.format(mcDir),    5,   'WWJets2L2Nu     {0}  2012'.format(mcTrigger)),
        cfg('ZZJets2L2Nu',     '{0}/ZZJetsTo2L2Nu'.format(mcDir),    5,   'ZZJets2L2Nu     {0}  2012'.format(mcTrigger)),
        cfg('ZZTo4e',          '{0}/ZZTo4e'.format(mcDir),           5,   'ZZ4e            {0}  2012'.format(mcTrigger)),
        cfg('ZZTo4mu',         '{0}/ZZTo4mu'.format(mcDir),          5,   'ZZ4mu           {0}  2012'.format(mcTrigger)),
        cfg('ZZTo4tau',        '{0}/ZZTo4tau'.format(mcDir),         5,   'ZZ4tau          {0}  2012'.format(mcTrigger)),
        cfg('ZZTo2e2mu',       '{0}/ZZTo2e2mu'.format(mcDir),        5,   'ZZ2e2mu         {0}  2012'.format(mcTrigger)),
        cfg('ZZTo2e2tau',      '{0}/ZZTo2e2tau'.format(mcDir),       5,   'ZZ2e2tau        {0}  2012'.format(mcTrigger)),
        cfg('ZZTo2mu2tau',     '{0}/ZZTo2mu2tau'.format(mcDir),      5,   'ZZ2mu2tau       {0}  2012'.format(mcTrigger))
    ])


inputSamples = []

if doData:
    inputSamples.extend(data)
if doBG:
    inputSamples.extend(bg)

if len(inputSamples) is not 0:
    batcher = b.BatchMaster(inputSamples, shortQueue = False, stageDir = 'batchStage', executable = executable, selection = selection + '_' + period)
    batcher.submit_to_batch()

