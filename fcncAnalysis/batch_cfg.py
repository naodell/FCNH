#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
lpchzgPath  = '/eos/uscms/store/user/lpchzg/'
filePath    = '/eos/uscms/store/user/naodell/data/'
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
        cfg('muon_2012A', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/DoubleMu_Run2012A', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/DoubleMu_Run2012B', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/DoubleMu_Run2012C', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/DoubleMu_Run2012D', 15, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012A', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012B', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012C', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/DoubleElectron_Run2012D', 15, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012A', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012B', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012C', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', lpchzgPath+'/nuTuples_v9.6_8TeV/Data/MuEG_Run2012D', 15, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        cfg('ZJets_M-50', filePath+'/nuTuples_v9_6_mc/DYJetsToLL_M-50', 30, 'ZJets_M-50 muon 2012'),
        cfg('ZJets_M-10To50', filePath+'/nuTuples_v9_6_mc/DYJetsToLL_M-10To50filter', 10, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZbbToLL', filePath+'/nuTuples_v9_6_mc/ZbbToLL', 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', filePath+'/nuTuples_v9_6_mc/ZGToLLG', 10, 'ZG muon 2012'),
        #cfg('WJets', filePath+'/nuTuples_v9_6_mc/WJetsToLNu', 20, 'WJets muon 2012'),
        #cfg('WGStarLNu2E', filePath+'/nuTuples_v9_6_mc/WGStarToLNu2E', 5, 'WGStarLNu2E muon 2012'),
        #cfg('WGStarLNu2Mu', filePath+'/nuTuples_v9_6_mc/WGStarToLNu2Mu', 5, 'WGStarLNu2Mu muon 2012'),
        #cfg('WGStarLNu2Tau', filePath+'/nuTuples_v9_6_mc/WGStarToLNu2Tau', 5, 'WGStarLNu2Tau muon 2012'),
        #cfg('WbbToLNu', filePath+'/nuTuples_v9_6_mc/WbbJetsToLNu', 20, 'WbbToLNu muon 2012'),
        #cfg('WbbToLNu', filePath+'/nuTuples_v9_6_mc/WbbToLL', 20, 'WbbToLNu muon 2012'),
        #cfg('WG', filePath+'/nuTuples_v9_6_mc/WGToLNuG', 10, 'WG muon 2012'),

        cfg('ttbarLep', filePath+'/nuTuples_v9_6_mc/TTJets', 30, 'ttbarLep muon 2012'),
        cfg('ttbarHad', filePath+'/nuTuples_v9_6_mc/TTJets', 30, 'ttbarHad muon 2012'),
        cfg('tbarW', filePath+'/nuTuples_v9_6_mc/Tbar_tW', 5, 'tbarW muon 2012'),
        cfg('tW', filePath+'/nuTuples_v9_6_mc/T_tW', 5, 'tW muon 2012'),
        cfg('t_t-channel', filePath+'/nuTuples_v9_6_mc/T_t', 5, 't_t-channel muon 2012'),
        cfg('tbar_t-channel', filePath+'/nuTuples_v9_6_mc/Tbar_t', 5, 'tbar_t-channel muon 2012'),
        cfg('ttW', filePath+'/nuTuples_v9_6_mc/TTWJets', 5, 'ttW muon 2012'),
        cfg('ttZ', filePath+'/nuTuples_v9_6_mc/TTZJets', 5, 'ttZ muon 2012'),
        cfg('ttG', filePath+'/nuTuples_v9_6_mc/TTGJets', 5, 'ttG muon 2012'),

        cfg('WWW', filePath+'/nuTuples_v9_6_mc/WWWJets', 5, 'WWW muon 2012'),
        cfg('WWZ', filePath+'/nuTuples_v9_6_mc/WWZNoGstarJets', 5, 'WWZ muon 2012'),
        cfg('WZZ', filePath+'/nuTuples_v9_6_mc/WZZNoGstarJets', 5, 'WZZ muon 2012'),
        cfg('ZZZ', filePath+'/nuTuples_v9_6_mc/ZZZNoGstarJets', 5, 'ZZZ muon 2012'),
        cfg('WWG', filePath+'/nuTuples_v9_6_mc/WWGJets', 5, 'WWG muon 2012'),

        cfg('ZZJets2L2Nu', filePath+'/nuTuples_v9_6_mc/ZZJetsTo2L2Nu', 5, 'ZZJets2L2Nu muon 2012'),
        #cfg('ZZJets2L2Q', filePath+'/nuTuples_v9_6_mc/ZZJetsTo2L2Q', 5, 'ZZJets2L2Q muon 2012'),
        cfg('ZZTo4e', filePath+'/nuTuples_v9_6_mc/ZZTo4e', 5, 'ZZ4e muon 2012'),
        cfg('ZZTo4mu', filePath+'/nuTuples_v9_6_mc/ZZTo4mu', 5, 'ZZ4mu muon 2012'),
        cfg('ZZTo4tau', filePath+'/nuTuples_v9_6_mc/ZZTo4tau', 5, 'ZZ4tau muon 2012'),
        cfg('ZZTo2e2mu', filePath+'/nuTuples_v9_6_mc/ZZTo2e2mu', 5, 'ZZ2e2mu muon 2012'),
        cfg('ZZTo2e2tau', filePath+'/nuTuples_v9_6_mc/ZZTo2e2tau', 5, 'ZZ2e2tau muon 2012'),
        cfg('ZZTo2mu2tau', filePath+'/nuTuples_v9_6_mc/ZZTo2mu2tau', 5, 'ZZ2mu2tau muon 2012'),
        cfg('WWJets2L2Nu', filePath+'/nuTuples_v9_6_mc/WWJetsTo2L2Nu', 5, 'WWJets2L2Nu muon 2012'),
        cfg('WZJets3LNu', filePath+'/nuTuples_v9_6_mc/WZJetsTo3LNu', 5, 'WZJets3LNu muon 2012'),
        cfg('WZJets2L2Q', filePath+'/nuTuples_v9_6_mc/WZJetsTo2L2Q', 10, 'WZJets2L2Q muon 2012'),

        #cfg('QCD_20-30_EM', filePath+'/nuTuples_v9_6_mc/QCD_Pt_20_30_EMEnriched', 10, 'QCD_20-30_EM muon 2012'),
        #cfg('QCD_30-80_EM', filePath+'/nuTuples_v9_6_mc/QCD_Pt_30_80_EMEnriched', 10, 'QCD_30-80_EM muon 2012'),
        #cfg('QCD_80-170_EM', filePath+'/nuTuples_v9_6_mc/QCD_Pt_80_170_EMEnriched', 10, 'QCD_80-170_EM muon 2012'),
        #cfg('QCD_170-250_EM', filePath+'/nuTuples_v9_6_mc/QCD_Pt_170_250_EMEnriched', 10, 'QCD_170-250_EM muon 2012'),
        #cfg('QCD_250-350_EM', filePath+'/nuTuples_v9_6_mc/QCD_Pt_250_350_EMEnriched', 10, 'QCD_250-350_EM muon 2012'),
        #cfg('QCD_350_EM', filePath+'/nuTuples_v9_6_mc/QCD_Pt_350_EMEnriched', 10, 'QCD_350_EM muon 2012'),
        #cfg('QCD_20_MU', filePath+'/nuTuples_v9_6_mc/QCD_Pt_20_MuEnrichedPt_15', 10, 'QCD_20_MU muon 2012'),

        cfg('ggHToZZ4L_M-125', filePath+'/nuTuples_v9_6_mc/GluGluToHToZZTo4L_M-125', 2, 'ggHToZZ4L_M-125 muon 2012'),
        cfg('ggHToWW2L2Nu_M-125', filePath+'/nuTuples_v9_6_mc/GluGluToHToWWTo2LAndTau2Nu_M-125', 2, 'ggHToWW2L2Nu_M-125 muon 2012'),
        cfg('WHToWWW3L_M-125', filePath+'/nuTuples_v9_6_mc/WH_HToWW_3l_M-125', 2, 'WHToWWW3L_M-125 muon 2012'),
        #cfg('TTH_M-125', filePath+'/nuTuples_v9_6_mc/TTH_Inclusive_M-125', 2, 'TTH_M-125 muon 2012')

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

