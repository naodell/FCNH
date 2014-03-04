#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
path  = '/eos/uscms/store/user/naodell/data/nuTuples_v7_4'
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
        cfg('muon_2012A', '{0}/DoubleMu_Run2012A'.format(path), 10, 'DATA_MUON muon 2012'),
        #cfg('muon_2012B', '{0}/DoubleMu_Run2012B'.format(path), 10, 'DATA_MUON muon 2012'),
        #cfg('muon_2012C', '{0}/DoubleMu_Run2012C'.format(path), 10, 'DATA_MUON muon 2012'),
        #cfg('muon_2012D', '{0}/DoubleMu_Run2012D'.format(path), 15, 'DATA_MUON muon 2012'),

        #cfg('electron_2012A', '{0}/DoubleElectron_Run2012A'.format(path), 10, 'DATA_ELECTRON electron 2012'),
        #cfg('electron_2012B', '{0}/DoubleElectron_Run2012B'.format(path), 10, 'DATA_ELECTRON electron 2012'),
        #cfg('electron_2012C', '{0}/DoubleElectron_Run2012C'.format(path), 10, 'DATA_ELECTRON electron 2012'),
        #cfg('electron_2012D', '{0}/DoubleElectron_Run2012D'.format(path), 15, 'DATA_ELECTRON electron 2012'),

        #cfg('muEG_2012A', '{0}/MuEG_Run2012A'.format(path), 10, 'DATA_MUEG muEG 2012'),
        #cfg('muEG_2012B', '{0}/MuEG_Run2012B'.format(path), 10, 'DATA_MUEG muEG 2012'),
        #cfg('muEG_2012C', '{0}/MuEG_Run2012C'.format(path), 10, 'DATA_MUEG muEG 2012'),
        #cfg('muEG_2012D', '{0}/MuEG_Run2012D'.format(path), 15, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        #cfg('ZJets_M-50', '{0}/DYJets'.format(path), 30, 'ZJets_M-50 muon 2012'),
        #cfg('ZJets_M-10To50', '{0}/DYJets_M-10To50'.format(path), 10, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZbbToLL', '{0}/ZbbToLL'.format(path), 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', '{0}/ZGToLLG'.format(path), 10, 'ZG muon 2012'),
        #cfg('WJets', '{0}/WJetsToLNu'.format(path), 20, 'WJets muon 2012'),
        #cfg('WGStarLNu2E', '{0}/WGStarToLNu2E'.format(path), 5, 'WGStarLNu2E muon 2012'),
        #cfg('WGStarLNu2Mu', '{0}/WGStarToLNu2Mu'.format(path), 5, 'WGStarLNu2Mu muon 2012'),
        #cfg('WGStarLNu2Tau', '{0}/WGStarToLNu2Tau'.format(path), 5, 'WGStarLNu2Tau muon 2012'),
        #cfg('WbbToLNu', '{0}/WbbJetsToLNu'.format(path), 20, 'WbbToLNu muon 2012'),
        #cfg('WbbToLNu', '{0}/WbbToLL'.format(path), 20, 'WbbToLNu muon 2012'),
        #cfg('WG', '{0}/WGToLNuG'.format(path), 10, 'WG muon 2012'),

        #cfg('ttbar', '{0}/TTJets'.format(path), 30, 'ttbar muon 2012'),
        #cfg('tbarW', '{0}/Tbar_tW'.format(path), 5, 'tbarW muon 2012'),
        #cfg('tW', '{0}/T_tW'.format(path), 5, 'tW muon 2012'),
        #cfg('t_t-channel', '{0}/T_t'.format(path), 5, 't_t-channel muon 2012'),
        #cfg('tbar_t-channel', '{0}/Tbar_t'.format(path), 5, 'tbar_t-channel muon 2012'),
        #cfg('ttW', '{0}/TTWJets'.format(path), 5, 'ttW muon 2012'),
        #cfg('ttZ', '{0}/TTZJets'.format(path), 5, 'ttZ muon 2012'),
        #cfg('ttG', '{0}/TTGJets'.format(path), 5, 'ttG muon 2012'),

        #cfg('WWW', '{0}/WWWJets'.format(path), 5, 'WWW muon 2012'),
        #cfg('WWZ', '{0}/WWZNoGstarJets'.format(path), 5, 'WWZ muon 2012'),
        #cfg('WZZ', '{0}/WZZNoGstarJets'.format(path), 5, 'WZZ muon 2012'),
        #cfg('ZZZ', '{0}/ZZZNoGstarJets'.format(path), 5, 'ZZZ muon 2012'),
        #cfg('WWG', '{0}/WWGJets'.format(path), 5, 'WWG muon 2012'),

        #cfg('ZZJets2L2Nu', '{0}/ZZJetsTo2L2Nu'.format(path), 5, 'ZZJets2L2Nu muon 2012'),
        #cfg('ZZJets2L2Q', '{0}/ZZJetsTo2L2Q'.format(path), 5, 'ZZJets2L2Q muon 2012'),
        #cfg('ZZJets4L', '{0}/ZZJetsTo4L'.format(path), 5, 'ZZJets4L muon 2012'),
        #cfg('ZZTo4e', '{0}/ZZTo4e'.format(path), 5, 'ZZ4e muon 2012'),
        #cfg('ZZTo4mu', '{0}/ZZTo4mu'.format(path), 5, 'ZZ4mu muon 2012'),
        #cfg('ZZTo4tau', '{0}/ZZTo4tau'.format(path), 5, 'ZZ4tau muon 2012'),
        #cfg('ZZTo2e2mu', '{0}/ZZTo2e2mu'.format(path), 5, 'ZZ2e2mu muon 2012'),
        #cfg('ZZTo2e2tau', '{0}/ZZTo2e2tau'.format(path), 5, 'ZZ2e2tau muon 2012'),
        #cfg('ZZTo2mu2tau', '{0}/ZZTo2mu2tau'.format(path), 5, 'ZZ2mu2tau muon 2012'),
        #cfg('WWJets2L2Nu', '{0}/WWJetsTo2L2Nu'.format(path), 5, 'WWJets2L2Nu muon 2012'),
        cfg('WZJets3LNu', '{0}/WZJetsTo3LNu'.format(path), 5, 'WZJets3LNu muon 2012'),
        #cfg('WZJets2L2Q', '{0}/WZJetsTo2L2Q'.format(path), 10, 'WZJets2L2Q muon 2012'),

        #cfg('QCD_20-30_EM', '{0}/QCD_Pt_20_30_EMEnriched'.format(path), 10, 'QCD_20-30_EM muon 2012'),
        #cfg('QCD_30-80_EM', '{0}/QCD_Pt_30_80_EMEnriched'.format(path), 10, 'QCD_30-80_EM muon 2012'),
        #cfg('QCD_80-170_EM', '{0}/QCD_Pt_80_170_EMEnriched'.format(path), 10, 'QCD_80-170_EM muon 2012'),
        #cfg('QCD_170-250_EM', '{0}/QCD_Pt_170_250_EMEnriched'.format(path), 10, 'QCD_170-250_EM muon 2012'),
        #cfg('QCD_250-350_EM', '{0}/QCD_Pt_250_350_EMEnriched'.format(path), 10, 'QCD_250-350_EM muon 2012'),
        #cfg('QCD_350_EM', '{0}/QCD_Pt_350_EMEnriched', 10, 'QCD_350_EM muon 2012'),
        #cfg('QCD_20_MU', '{0}/QCD_Pt_20_MuEnrichedPt_15', 10, 'QCD_20_MU muon 2012'),

        #cfg('ggHToZZ4L_M-125', '{0}/GluGluToHToZZTo4L_M-125', 2, 'ggHToZZ4L_M-125 muon 2012'),
        #cfg('ggHToWW2L2Nu_M-125', '{0}/GluGluToHToWWTo2LAndTau2Nu_M-125', 2, 'ggHToWW2L2Nu_M-125 muon 2012'),
        #cfg('WHToWWW3L_M-125', '{0}/WH_HToWW_3l_M-125', 2, 'WHToWWW3L_M-125 muon 2012'),
        #cfg('TTH_M-125', '{0}/TTH_Inclusive_M-125', 2, 'TTH_M-125 muon 2012')
        ])

    signal.extend([
        #cfg('FCNC_M125_tHj', filePath+'/nuTuples_v7_4/TToFCNHToWWTo2l2nuPlusTop_M125', 1, 'FCNC_M125_t mc 2012'),
        #cfg('FCNC_M125_tbarHj', filePath+'/nuTuples_v7_4/TbarToFCNHToWWTo2l2nuPlusTop_M125/', 1, 'FCNC_M125_tbar mc 2012')
        cfg('FCNC_M125_tHj', filePath+'/nuTuples_v7_4/T_FCNH_M-125_WW_lep', 1, 'FCNC_M125_t mc 2012'),
        #cfg('FCNC_M125_tHj', filePath+'/nuTuples_v6_8TeV/FCNH_M125_t', 1, 'FCNC_M125_t mc 2012'),
        #cfg('FCNC_M125_tbarHj', filePath+'/nuTuples_v6_8TeV/FCNH_M125_tbar', 1, 'FCNC_M125_tbar mc 2012'),
        cfg('FCNC_M125_tHj_semilep', filePath+'/nuTuples_v7_4/T_FCNH_M-125_WW_semihadronic', 1, 'FCNC_M125_t_semilep mc 2012'),
        cfg('FCNC_M125_tHj_ZZ', filePath+'/nuTuples_v7_4/T_FCNH_M-125_ZZ', 1, 'FCNC_M125_t_ZZ mc 2012'),
        cfg('FCNC_M125_tHj_TauTau', filePath+'/nuTuples_v7_4/T_FCNH_M-125_TauTau', 1, 'FCNC_M125_t_TauTau mc 2012')
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

