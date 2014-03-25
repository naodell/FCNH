#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
mcDir       = '/tthome/naodell/storage/data/nuTuples_v7_4'
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
        cfg('muon_2012A', '{0}/DoubleMu_Run2012A'.format(path), 20, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', '{0}/DoubleMu_Run2012B'.format(path), 20, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', '{0}/DoubleMu_Run2012C'.format(path), 20, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', '{0}/DoubleMu_Run2012D'.format(path), 25, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', '{0}/DoubleElectron_Run2012A'.format(path), 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', '{0}/DoubleElectron_Run2012B'.format(path), 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', '{0}/DoubleElectron_Run2012C'.format(path), 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', '{0}/DoubleElectron_Run2012D'.format(path), 25, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', '{0}/MuEG_Run2012A'.format(path), 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', '{0}/MuEG_Run2012B'.format(path), 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', '{0}/MuEG_Run2012C'.format(path), 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', '{0}/MuEG_Run2012D'.format(path), 25, 'DATA_MUEG muEG 2012')
    ])

    bg.extend([
        cfg('ZJets_M-50',      '{0}/DYJets'.format(path),           30,  'ZJets_M-50      muon  2012'),
        cfg('ZJets_M-10To50',  '{0}/DYJets_M-10To50'.format(path),  10,  'ZJets_M-10To50  muon  2012'),
        cfg('WbbToLNu',        '{0}/WbbToLL'.format(path),          20,  'WbbToLNu        muon  2012'),
        cfg('WjetToLNu',       '{0}/WJetsToLNu'.format(path),       20,  'WJetsToLNu      muon  2012'),

        cfg('ttbar',           '{0}/TTJets'.format(path),           30,  'ttbar           muon  2012'),
        cfg('tbarW',           '{0}/Tbar_tW'.format(path),          5,   'tbarW           muon  2012'),
        cfg('tW',              '{0}/T_tW'.format(path),             5,   'tW              muon  2012'),

        cfg('WZJets3LNu',      '{0}/WZJetsTo3LNu'.format(path),     5,   'WZJets3LNu      muon  2012'),
        cfg('WWJets2L2Nu',     '{0}/WWJetsTo2L2Nu'.format(path),    5,   'WWJets2L2Nu     muon  2012'),
        cfg('ZZJets2L2Nu',     '{0}/ZZJetsTo2L2Nu'.format(path),    5,   'ZZJets2L2Nu     muon  2012'),
        cfg('ZZTo4e',          '{0}/ZZTo4e'.format(path),           5,   'ZZ4e            muon  2012'),
        cfg('ZZTo4mu',         '{0}/ZZTo4mu'.format(path),          5,   'ZZ4mu           muon  2012'),
        cfg('ZZTo4tau',        '{0}/ZZTo4tau'.format(path),         5,   'ZZ4tau          muon  2012'),
        cfg('ZZTo2e2mu',       '{0}/ZZTo2e2mu'.format(path),        5,   'ZZ2e2mu         muon  2012'),
        cfg('ZZTo2e2tau',      '{0}/ZZTo2e2tau'.format(path),       5,   'ZZ2e2tau        muon  2012'),
        cfg('ZZTo2mu2tau',     '{0}/ZZTo2mu2tau'.format(path),      5,   'ZZ2mu2tau       muon  2012')
    ])


inputSamples = []

if doData:
    inputSamples.extend(data)
if doBG:
    inputSamples.extend(bg)

if len(inputSamples) is not 0:
    batcher = b.BatchMaster(inputSamples, shortQueue = False, stageDir = 'batchStage', executable = executable, selection = selection + '_' + period)
    batcher.submit_to_batch()

