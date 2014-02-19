#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
#inputDir      = '/tthome/naodell/storage/data/nuTuples_v7_4'
inputDir      = '/tthome/naodell/storage/data/nuTuples_v9.6_8TeV/MC'
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
        cfg('muon_2012A', inputDir+'/DoubleMu_Run2012A', 40, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', inputDir+'/DoubleMu_Run2012B', 40, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', inputDir+'/DoubleMu_Run2012C', 40, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', inputDir+'/DoubleMu_Run2012D', 45, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', inputDir+'/DoubleElectron_Run2012A', 40, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', inputDir+'/DoubleElectron_Run2012B', 40, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', inputDir+'/DoubleElectron_Run2012C', 40, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', inputDir+'/DoubleElectron_Run2012D', 45, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', inputDir+'/MuEG_Run2012A', 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', inputDir+'/MuEG_Run2012B', 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', inputDir+'/MuEG_Run2012C', 20, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', inputDir+'/MuEG_Run2012D', 25, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        cfg('ZJets_M-50', inputDir+'/DYJetsToLL_M-50', 50, 'ZJets_M-50 muon 2012'),
        cfg('ZJets_M-10To50', inputDir+'/DYJetsToLL_M-10To50filter', 35, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZJets_M-50', inputDir+'/DYJets', 50, 'ZJets_M-50 muon 2012'),
        #cfg('ZJets_M-10To50', inputDir+'/DYJets_M-10To50', 35, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZbbToLL', inputDir+'/ZbbToLL', 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', inputDir+'/ZGToLLG', 10, 'ZG muon 2012'),
        #cfg('WJets', inputDir+'/WJetsToLNu', 20, 'WJets muon 2012'),
        #cfg('WGStarLNu2E', inputDir+'/WGStarToLNu2E', 5, 'WGStarLNu2E muon 2012'),
        #cfg('WGStarLNu2Mu', inputDir+'/WGStarToLNu2Mu', 5, 'WGStarLNu2Mu muon 2012'),
        #cfg('WGStarLNu2Tau', inputDir+'/WGStarToLNu2Tau', 5, 'WGStarLNu2Tau muon 2012'),
        #cfg('WbbToLNu', inputDir+'/WbbJetsToLNu', 20, 'WbbToLNu muon 2012'),
        #cfg('WbbToLNu', inputDir+'/WbbToLL', 20, 'WbbToLNu muon 2012'),
        #cfg('WG', inputDir+'/WGToLNuG', 10, 'WG muon 2012'),

        #cfg('ttbar', inputDir+'/TTJets', 30, 'ttbar muon 2012'),
        #cfg('tbarW', inputDir+'/Tbar_tW', 5, 'tbarW muon 2012'),
        #cfg('tW', inputDir+'/T_tW', 5, 'tW muon 2012'),
        #cfg('t_t-channel', inputDir+'/T_t', 5, 't_t-channel muon 2012'),
        #cfg('tbar_t-channel', inputDir+'/Tbar_t', 5, 'tbar_t-channel muon 2012'),
        #cfg('ttW', inputDir+'/TTWJets', 5, 'ttW muon 2012'),
        #cfg('ttZ', inputDir+'/TTZJets', 5, 'ttZ muon 2012'),
        #cfg('ttG', inputDir+'/TTGJets', 5, 'ttG muon 2012'),

        #cfg('WWW', inputDir+'/WWWJets', 5, 'WWW muon 2012'),
        #cfg('WWZ', inputDir+'/WWZNoGstarJets', 5, 'WWZ muon 2012'),
        #cfg('WZZ', inputDir+'/WZZNoGstarJets', 5, 'WZZ muon 2012'),
        #cfg('ZZZ', inputDir+'/ZZZNoGstarJets', 5, 'ZZZ muon 2012'),
        #cfg('WWG', inputDir+'/WWGJets', 5, 'WWG muon 2012'),

        #cfg('ZZJets2L2Nu', inputDir+'/ZZJetsTo2L2Nu', 5, 'ZZJets2L2Nu muon 2012'),
        #cfg('ZZJets2L2Q', inputDir+'/ZZJetsTo2L2Q', 5, 'ZZJets2L2Q muon 2012'),
        #cfg('ZZJets4L', inputDir+'/ZZJetsTo4L', 5, 'ZZJets4L muon 2012'),
        #cfg('ZZTo4e', inputDir+'/ZZTo4e', 5, 'ZZ4e muon 2012'),
        #cfg('ZZTo4mu', inputDir+'/ZZTo4mu', 5, 'ZZ4mu muon 2012'),
        #cfg('ZZTo4tau', inputDir+'/ZZTo4tau', 5, 'ZZ4tau muon 2012'),
        #cfg('ZZTo2e2mu', inputDir+'/ZZTo2e2mu', 5, 'ZZ2e2mu muon 2012'),
        #cfg('ZZTo2e2tau', inputDir+'/ZZTo2e2tau', 5, 'ZZ2e2tau muon 2012'),
        #cfg('ZZTo2mu2tau', inputDir+'/ZZTo2mu2tau', 5, 'ZZ2mu2tau muon 2012'),
        #cfg('WWJets2L2Nu', inputDir+'/WWJetsTo2L2Nu', 5, 'WWJets2L2Nu muon 2012'),
        #cfg('WZJets3LNu', inputDir+'/WZJetsTo3LNu', 5, 'WZJets3LNu muon 2012'),
        #cfg('WZJets2L2Q', inputDir+'/WZJetsTo2L2Q', 10, 'WZJets2L2Q muon 2012'),

        #cfg('QCD_20-30_EM', inputDir+'/QCD_Pt_20_30_EMEnriched', 10, 'QCD_20-30_EM muon 2012'),
        #cfg('QCD_30-80_EM', inputDir+'/QCD_Pt_30_80_EMEnriched', 10, 'QCD_30-80_EM muon 2012'),
        #cfg('QCD_80-170_EM', inputDir+'/QCD_Pt_80_170_EMEnriched', 10, 'QCD_80-170_EM muon 2012'),
        #cfg('QCD_170-250_EM', inputDir+'/QCD_Pt_170_250_EMEnriched', 10, 'QCD_170-250_EM muon 2012'),
        #cfg('QCD_250-350_EM', inputDir+'/QCD_Pt_250_350_EMEnriched', 10, 'QCD_250-350_EM muon 2012'),
        #cfg('QCD_350_EM', inputDir+'/QCD_Pt_350_EMEnriched', 10, 'QCD_350_EM muon 2012'),
        #cfg('QCD_20_MU', inputDir+'/QCD_Pt_20_MuEnrichedPt_15', 10, 'QCD_20_MU muon 2012'),

        #cfg('ggHToZZ4L_M-125', inputDir+'/GluGluToHToZZTo4L_M-125', 2, 'ggHToZZ4L_M-125 muon 2012'),
        #cfg('ggHToWW2L2Nu_M-125', inputDir+'/GluGluToHToWWTo2LAndTau2Nu_M-125', 2, 'ggHToWW2L2Nu_M-125 muon 2012'),
        #cfg('WHToWWW3L_M-125', inputDir+'/WH_HToWW_3l_M-125', 2, 'WHToWWW3L_M-125 muon 2012'),
        #cfg('TTH_M-125', inputDir+'/TTH_Inclusive_M-125', 2, 'TTH_M-125 muon 2012')

        ])

    signal.extend([
        #cfg('FCNC_M125_tHj', inputDir+'/nuTuples_v7_4/TToFCNHToWWTo2l2nuPlusTop_M125', 1, 'FCNC_M125_t mc 2012'),
        #cfg('FCNC_M125_tbarHj', inputDir+'/nuTuples_v7_4/TbarToFCNHToWWTo2l2nuPlusTop_M125/', 1, 'FCNC_M125_tbar mc 2012')
        cfg('FCNC_M125_tHj', inputDir+'/nuTuples_v6_8TeV/FCNH_M125_t', 1, 'FCNC_M125_t mc 2012'),
        cfg('FCNC_M125_tbarHj', inputDir+'/nuTuples_v6_8TeV/FCNH_M125_tbar', 1, 'FCNC_M125_tbar mc 2012')
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

