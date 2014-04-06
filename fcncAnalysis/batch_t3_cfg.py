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
        cfg('muon_2012A', '{0}/DoubleMu_Run2012A'.format(dataDir), 50, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', '{0}/DoubleMu_Run2012B'.format(dataDir), 50, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', '{0}/DoubleMu_Run2012C'.format(dataDir), 50, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', '{0}/DoubleMu_Run2012D'.format(dataDir), 55, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', '{0}/DoubleElectron_Run2012A'.format(dataDir), 50, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', '{0}/DoubleElectron_Run2012B'.format(dataDir), 50, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', '{0}/DoubleElectron_Run2012C'.format(dataDir), 50, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', '{0}/DoubleElectron_Run2012D'.format(dataDir), 55, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', '{0}/MuEG_Run2012A'.format(dataDir), 30, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', '{0}/MuEG_Run2012B'.format(dataDir), 30, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', '{0}/MuEG_Run2012C'.format(dataDir), 30, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', '{0}/MuEG_Run2012D'.format(dataDir), 35, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        cfg('ZJets_M-50', mcDir+'/DYJetsToLL_M-50', 50, 'ZJets_M-50 muon 2012'),
        cfg('ZJets_M-10To50', mcDir+'/DYJetsToLL_M-10To50filter', 35, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZbbToLL', mcDir+'/ZbbToLL', 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', mcDir+'/ZGToLLG', 10, 'ZG muon 2012'),
        #cfg('WJets', mcDir+'/WJetsToLNu', 20, 'WJets muon 2012'),
        #cfg('WGStarLNu2E', mcDir+'/WGstarToLNu2E', 5, 'WGStarLNu2E muon 2012'),
        #cfg('WGStarLNu2Mu', mcDir+'/WGstarToLNu2Mu', 5, 'WGStarLNu2Mu muon 2012'),
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
        cfg('ttW', mcDir+'/TTWJets', 5, 'ttW muon 2012'),
        cfg('ttZ', mcDir+'/TTZJets', 5, 'ttZ muon 2012'),
        cfg('ttG', mcDir+'/TTGJets', 5, 'ttG muon 2012'),

        #cfg('WWW', mcDir+'/WWWJets', 5, 'WWW muon 2012'),
        #cfg('WWZ', mcDir+'/WWZNoGstarJets', 5, 'WWZ muon 2012'),
        #cfg('WZZ', mcDir+'/WZZNoGstarJets', 5, 'WZZ muon 2012'),
        #cfg('ZZZ', mcDir+'/ZZZNoGstarJets', 5, 'ZZZ muon 2012'),
        #cfg('WWG', mcDir+'/WWGJets', 5, 'WWG muon 2012'),

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

        #cfg('QCD_20-30_EM', mcDir+'/QCD_Pt_20_30_EMEnriched', 10, 'QCD_20-30_EM muon 2012'),
        #cfg('QCD_30-80_EM', mcDir+'/QCD_Pt_30_80_EMEnriched', 10, 'QCD_30-80_EM muon 2012'),
        #cfg('QCD_80-170_EM', mcDir+'/QCD_Pt_80_170_EMEnriched', 10, 'QCD_80-170_EM muon 2012'),
        #cfg('QCD_170-250_EM', mcDir+'/QCD_Pt_170_250_EMEnriched', 10, 'QCD_170-250_EM muon 2012'),
        #cfg('QCD_250-350_EM', mcDir+'/QCD_Pt_250_350_EMEnriched', 10, 'QCD_250-350_EM muon 2012'),
        #cfg('QCD_350_EM', mcDir+'/QCD_Pt_350_EMEnriched', 10, 'QCD_350_EM muon 2012'),
        #cfg('QCD_20_MU', mcDir+'/QCD_Pt_20_MuEnrichedPt_15', 10, 'QCD_20_MU muon 2012'),

        #cfg('ggHToZZ4L_M-125', mcDir+'/GluGluToHToZZTo4L_M-125', 5, 'ggHToZZ4L_M-125 muon 2012'),
        #cfg('ggHToWW2L2Nu_M-125', mcDir+'/GluGluToHToWWTo2LAndTau2Nu_M-125', 5, 'ggHToWW2L2Nu_M-125 muon 2012'),
        #cfg('WHToWWW3L_M-125', mcDir+'/WH_HToWW_3l_M-125', 5, 'WHToWWW3L_M-125 muon 2012'),
        #cfg('TTH_M-125', mcDir+'/TTH_Inclusive_M-125', 5, 'TTH_M-125 muon 2012')

        ])

    signal.extend([
        #cfg('FCNC_M125_tHj', mcDir+'/TToFCNHToWWTo2l2nuPlusTop_M125', 1, 'FCNC_M125_t mc 2012'),
        #cfg('FCNC_M125_tbarHj', mcDir+'/TbarToFCNHToWWTo2l2nuPlusTop_M125/', 1, 'FCNC_M125_tbar mc 2012')
        cfg('FCNC_M125_tHj', mcDir+'/FCNH_M125_t', 1, 'FCNC_M125_t mc 2012'),
        cfg('FCNC_M125_tbarHj', mcDir+'/FCNH_M125_tbar', 1, 'FCNC_M125_tbar mc 2012')
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

