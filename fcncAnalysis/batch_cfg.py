#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
inputDir      = '/tthome/naodell/storage/data'
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
    #data.extend([
    #    cfg('muon_2012A', dCache+'/naodell/nuTuples_v7_4/DoubleMu_Run2012A', 10, 'DATA_MUON muon 2012'),
    #    cfg('muon_2012B', dCache+'/naodell/nuTuples_v7_4/DoubleMu_Run2012B', 10, 'DATA_MUON muon 2012'),
    #    cfg('muon_2012C', dCache+'/naodell/nuTuples_v7_4/DoubleMu_Run2012C', 10, 'DATA_MUON muon 2012'),
    #    cfg('muon_2012D', dCache+'/naodell/nuTuples_v7_4/DoubleMu_Run2012D', 15, 'DATA_MUON muon 2012'),

    #    cfg('electron_2012A', dCache+'/naodell/nuTuples_v7_4/DoubleElectron_Run2012A', 10, 'DATA_ELECTRON electron 2012'),
    #    cfg('electron_2012B', dCache+'/naodell/nuTuples_v7_4/DoubleElectron_Run2012B', 10, 'DATA_ELECTRON electron 2012'),
    #    cfg('electron_2012C', dCache+'/naodell/nuTuples_v7_4/DoubleElectron_Run2012C', 10, 'DATA_ELECTRON electron 2012'),
    #    cfg('electron_2012D', dCache+'/naodell/nuTuples_v7_4/DoubleElectron_Run2012D', 15, 'DATA_ELECTRON electron 2012'),

    #    cfg('muEG_2012A', dCache+'/naodell/nuTuples_v7_4/MuEG_Run2012A', 10, 'DATA_MUEG muEG 2012'),
    #    cfg('muEG_2012B', dCache+'/naodell/nuTuples_v7_4/MuEG_Run2012B', 10, 'DATA_MUEG muEG 2012'),
    #    cfg('muEG_2012C', dCache+'/naodell/nuTuples_v7_4/MuEG_Run2012C', 10, 'DATA_MUEG muEG 2012'),
    #    cfg('muEG_2012D', dCache+'/naodell/nuTuples_v7_4/MuEG_Run2012D', 15, 'DATA_MUEG muEG 2012')
    #    ])

    data.extend([
<<<<<<< HEAD
        cfg('muon_2012A', inputDir+'/nuTuples_v7_4/DoubleMu_Run2012A', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', inputDir+'/nuTuples_v7_4/DoubleMu_Run2012B', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', inputDir+'/nuTuples_v7_4/DoubleMu_Run2012C', 10, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', inputDir+'/nuTuples_v7_4/DoubleMu_Run2012D', 15, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', inputDir+'/nuTuples_v7_4/DoubleElectron_Run2012A', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', inputDir+'/nuTuples_v7_4/DoubleElectron_Run2012B', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', inputDir+'/nuTuples_v7_4/DoubleElectron_Run2012C', 10, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', inputDir+'/nuTuples_v7_4/DoubleElectron_Run2012D', 15, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', inputDir+'/nuTuples_v7_4/MuEG_Run2012A', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', inputDir+'/nuTuples_v7_4/MuEG_Run2012B', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', inputDir+'/nuTuples_v7_4/MuEG_Run2012C', 10, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', inputDir+'/nuTuples_v7_4/MuEG_Run2012D', 15, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        cfg('ZJets_M-50', inputDir+'/nuTuples_v7_4/DYJets', 30, 'ZJets_M-50 muon 2012')
        #cfg('ZJets_M-10To50', inputDir+'/nuTuples_v7_4/DYJets_M-10To50', 10, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZbbToLL', inputDir+'/nuTuples_v7_4/ZbbToLL', 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', inputDir+'/nuTuples_v7_4/ZGToLLG', 10, 'ZG muon 2012'),
        #cfg('WJets', inputDir+'/nuTuples_v7_4/WJetsToLNu', 20, 'WJets muon 2012'),
        #cfg('WGStarLNu2E', inputDir+'/nuTuples_v7_4/WGStarToLNu2E', 5, 'WGStarLNu2E muon 2012'),
        #cfg('WGStarLNu2Mu', inputDir+'/nuTuples_v7_4/WGStarToLNu2Mu', 5, 'WGStarLNu2Mu muon 2012'),
        #cfg('WGStarLNu2Tau', inputDir+'/nuTuples_v7_4/WGStarToLNu2Tau', 5, 'WGStarLNu2Tau muon 2012'),
        #cfg('WbbToLNu', inputDir+'/nuTuples_v7_4/WbbJetsToLNu', 20, 'WbbToLNu muon 2012'),
        #cfg('WbbToLNu', inputDir+'/nuTuples_v7_4/WbbToLL', 20, 'WbbToLNu muon 2012'),
        #cfg('WG', inputDir+'/nuTuples_v7_4/WGToLNuG', 10, 'WG muon 2012'),

        #cfg('ttbar', inputDir+'/nuTuples_v7_4/TTJets', 30, 'ttbar muon 2012'),
        #cfg('tbarW', inputDir+'/nuTuples_v7_4/Tbar_tW', 5, 'tbarW muon 2012'),
        #cfg('tW', inputDir+'/nuTuples_v7_4/T_tW', 5, 'tW muon 2012'),
        #cfg('t_t-channel', inputDir+'/nuTuples_v7_4/T_t', 5, 't_t-channel muon 2012'),
        #cfg('tbar_t-channel', inputDir+'/nuTuples_v7_4/Tbar_t', 5, 'tbar_t-channel muon 2012'),
        #cfg('ttW', inputDir+'/nuTuples_v7_4/TTWJets', 5, 'ttW muon 2012'),
        #cfg('ttZ', inputDir+'/nuTuples_v7_4/TTZJets', 5, 'ttZ muon 2012'),
        #cfg('ttG', inputDir+'/nuTuples_v7_4/TTGJets', 5, 'ttG muon 2012'),

        #cfg('WWW', inputDir+'/nuTuples_v7_4/WWWJets', 5, 'WWW muon 2012'),
        #cfg('WWZ', inputDir+'/nuTuples_v7_4/WWZNoGstarJets', 5, 'WWZ muon 2012'),
        #cfg('WZZ', inputDir+'/nuTuples_v7_4/WZZNoGstarJets', 5, 'WZZ muon 2012'),
        #cfg('ZZZ', inputDir+'/nuTuples_v7_4/ZZZNoGstarJets', 5, 'ZZZ muon 2012'),
        #cfg('WWG', inputDir+'/nuTuples_v7_4/WWGJets', 5, 'WWG muon 2012'),

        #cfg('ZZJets2L2Nu', inputDir+'/nuTuples_v7_4/ZZJetsTo2L2Nu', 5, 'ZZJets2L2Nu muon 2012'),
        #cfg('ZZJets2L2Q', inputDir+'/nuTuples_v7_4/ZZJetsTo2L2Q', 5, 'ZZJets2L2Q muon 2012'),
        #cfg('ZZJets4L', inputDir+'/nuTuples_v7_4/ZZJetsTo4L', 5, 'ZZJets4L muon 2012'),
        #cfg('ZZTo4e', inputDir+'/nuTuples_v7_4/ZZTo4e', 5, 'ZZ4e muon 2012'),
        #cfg('ZZTo4mu', inputDir+'/nuTuples_v7_4/ZZTo4mu', 5, 'ZZ4mu muon 2012'),
        #cfg('ZZTo4tau', inputDir+'/nuTuples_v7_4/ZZTo4tau', 5, 'ZZ4tau muon 2012'),
        #cfg('ZZTo2e2mu', inputDir+'/nuTuples_v7_4/ZZTo2e2mu', 5, 'ZZ2e2mu muon 2012'),
        #cfg('ZZTo2e2tau', inputDir+'/nuTuples_v7_4/ZZTo2e2tau', 5, 'ZZ2e2tau muon 2012'),
        #cfg('ZZTo2mu2tau', inputDir+'/nuTuples_v7_4/ZZTo2mu2tau', 5, 'ZZ2mu2tau muon 2012'),
        #cfg('WWJets2L2Nu', inputDir+'/nuTuples_v7_4/WWJetsTo2L2Nu', 5, 'WWJets2L2Nu muon 2012'),
        #cfg('WZJets3LNu', inputDir+'/nuTuples_v7_4/WZJetsTo3LNu', 5, 'WZJets3LNu muon 2012'),
        #cfg('WZJets2L2Q', inputDir+'/nuTuples_v7_4/WZJetsTo2L2Q', 10, 'WZJets2L2Q muon 2012'),

        #cfg('QCD_20-30_EM', inputDir+'/nuTuples_v7_4/QCD_Pt_20_30_EMEnriched', 10, 'QCD_20-30_EM muon 2012'),
        #cfg('QCD_30-80_EM', inputDir+'/nuTuples_v7_4/QCD_Pt_30_80_EMEnriched', 10, 'QCD_30-80_EM muon 2012'),
        #cfg('QCD_80-170_EM', inputDir+'/nuTuples_v7_4/QCD_Pt_80_170_EMEnriched', 10, 'QCD_80-170_EM muon 2012'),
        #cfg('QCD_170-250_EM', inputDir+'/nuTuples_v7_4/QCD_Pt_170_250_EMEnriched', 10, 'QCD_170-250_EM muon 2012'),
        #cfg('QCD_250-350_EM', inputDir+'/nuTuples_v7_4/QCD_Pt_250_350_EMEnriched', 10, 'QCD_250-350_EM muon 2012'),
        #cfg('QCD_350_EM', inputDir+'/nuTuples_v7_4/QCD_Pt_350_EMEnriched', 10, 'QCD_350_EM muon 2012'),
        #cfg('QCD_20_MU', inputDir+'/nuTuples_v7_4/QCD_Pt_20_MuEnrichedPt_15', 10, 'QCD_20_MU muon 2012'),

        #cfg('ggHToZZ4L_M-125', inputDir+'/nuTuples_v7_4/GluGluToHToZZTo4L_M-125', 2, 'ggHToZZ4L_M-125 muon 2012'),
        #cfg('ggHToWW2L2Nu_M-125', inputDir+'/nuTuples_v7_4/GluGluToHToWWTo2LAndTau2Nu_M-125', 2, 'ggHToWW2L2Nu_M-125 muon 2012'),
        #cfg('WHToWWW3L_M-125', inputDir+'/nuTuples_v7_4/WH_HToWW_3l_M-125', 2, 'WHToWWW3L_M-125 muon 2012'),
        #cfg('TTH_M-125', inputDir+'/nuTuples_v7_4/TTH_Inclusive_M-125', 2, 'TTH_M-125 muon 2012')

        ])

    signal.extend([
        #cfg('FCNC_M125_tHj', inputDir+'/nuTuples_v7_4/TToFCNHToWWTo2l2nuPlusTop_M125', 1, 'FCNC_M125_t mc 2012'),
        #cfg('FCNC_M125_tbarHj', inputDir+'/nuTuples_v7_4/TbarToFCNHToWWTo2l2nuPlusTop_M125/', 1, 'FCNC_M125_tbar mc 2012')
        cfg('FCNC_M125_tHj', inputDir+'/nuTuples_v6_8TeV/FCNH_M125_t', 1, 'FCNC_M125_t mc 2012'),
        cfg('FCNC_M125_tbarHj', inputDir+'/nuTuples_v6_8TeV/FCNH_M125_tbar', 1, 'FCNC_M125_tbar mc 2012')
=======
        #cfg('muon_2012A', dCache+'/naodell/nuTuples_v6_8TeV/DoubleMu_Run2012A', 10, 'DATA_MUON muon 2012'),
        #cfg('muon_2012B', dCache+'/naodell/nuTuples_v6_8TeV/DoubleMu_Run2012B', 10, 'DATA_MUON muon 2012'),
        #cfg('muon_2012C', dCache+'/naodell/nuTuples_v6_8TeV/DoubleMu_Run2012C', 10, 'DATA_MUON muon 2012'),
        #cfg('muon_2012D', dCache+'/naodell/nuTuples_v6_8TeV/DoubleMu_Run2012D', 15, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', dCache+'/naodell/nuTuples_v6_8TeV/DoubleElectron_Run2012A', 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', dCache+'/naodell/nuTuples_v6_8TeV/DoubleElectron_Run2012B', 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', dCache+'/naodell/nuTuples_v6_8TeV/DoubleElectron_Run2012C', 20, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', dCache+'/naodell/nuTuples_v6_8TeV/DoubleElectron_Run2012D', 25, 'DATA_ELECTRON electron 2012')

        #cfg('muEG_2012A', dCache+'/naodell/nuTuples_v6_8TeV/MuEG_Run2012A', 10, 'DATA_MUEG muEG 2012'),
        #cfg('muEG_2012B', dCache+'/naodell/nuTuples_v6_8TeV/MuEG_Run2012B', 10, 'DATA_MUEG muEG 2012'),
        #cfg('muEG_2012C', dCache+'/naodell/nuTuples_v6_8TeV/MuEG_Run2012C', 10, 'DATA_MUEG muEG 2012'),
        #cfg('muEG_2012D', dCache+'/naodell/nuTuples_v6_8TeV/MuEG_Run2012D', 15, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        cfg('ZJets_M-50', dCache+'/naodell/nuTuples_v7_4/DYJets', 30, 'ZJets_M-50 muon 2012'),
        cfg('ZJets_M-10To50', dCache+'/naodell/nuTuples_v7_4/DYJets_M-10To50', 10, 'ZJets_M-10To50 muon 2012'),
        #cfg('ZbbToLL', dCache+'/naodell/nuTuples_v7_4/ZbbToLL', 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', dCache+'/naodell/nuTuples_v7_4/ZGToLLG', 10, 'ZG muon 2012'),
        #cfg('WJets', dCache+'/naodell/nuTuples_v7_4/WJetsToLNu', 20, 'WJets muon 2012'),
        #cfg('WGStarLNu2E', dCache+'/naodell/nuTuples_v7_4/WGStarToLNu2E', 5, 'WGStarLNu2E muon 2012'),
        #cfg('WGStarLNu2Mu', dCache+'/naodell/nuTuples_v7_4/WGStarToLNu2Mu', 5, 'WGStarLNu2Mu muon 2012'),
        #cfg('WGStarLNu2Tau', dCache+'/naodell/nuTuples_v7_4/WGStarToLNu2Tau', 5, 'WGStarLNu2Tau muon 2012'),
        #cfg('WbbToLNu', dCache+'/naodell/nuTuples_v7_4/WbbJetsToLNu', 20, 'WbbToLNu muon 2012'),
<<<<<<< HEAD
=======
        #cfg('WbbToLNu', dCache+'/naodell/nuTuples_v7_4/WbbToLL', 20, 'WbbToLNu muon 2012'),
>>>>>>> recover
        #cfg('WG', dCache+'/naodell/nuTuples_v7_4/WGToLNuG', 10, 'WG muon 2012'),

        cfg('ttbar', dCache+'/naodell/nuTuples_v7_4/TTJets', 30, 'ttbar muon 2012'),
        cfg('tbarW', dCache+'/naodell/nuTuples_v7_4/Tbar_tW', 5, 'tbarW muon 2012'),
        cfg('tW', dCache+'/naodell/nuTuples_v7_4/T_tW', 5, 'tW muon 2012'),
        cfg('t_t-channel', dCache+'/naodell/nuTuples_v7_4/T_t', 5, 't_t-channel muon 2012'),
        cfg('tbar_t-channel', dCache+'/naodell/nuTuples_v7_4/Tbar_t', 5, 'tbar_t-channel muon 2012'),
        cfg('ttW', dCache+'/naodell/nuTuples_v7_4/TTWJets', 5, 'ttW muon 2012'),
        cfg('ttZ', dCache+'/naodell/nuTuples_v7_4/TTZJets', 5, 'ttZ muon 2012'),
        cfg('ttG', dCache+'/naodell/nuTuples_v7_4/TTGJets', 5, 'ttG muon 2012'),

        cfg('WWW', dCache+'/naodell/nuTuples_v7_4/WWWJets', 5, 'WWW muon 2012'),
        cfg('WWZ', dCache+'/naodell/nuTuples_v7_4/WWZNoGstarJets', 5, 'WWZ muon 2012'),
        cfg('WZZ', dCache+'/naodell/nuTuples_v7_4/WZZNoGstarJets', 5, 'WZZ muon 2012'),
        cfg('ZZZ', dCache+'/naodell/nuTuples_v7_4/ZZZNoGstarJets', 5, 'ZZZ muon 2012'),
        cfg('WWG', dCache+'/naodell/nuTuples_v7_4/WWGJets', 5, 'WWG muon 2012'),

        cfg('ZZJets2L2Nu', dCache+'/naodell/nuTuples_v7_4/ZZJetsTo2L2Nu', 5, 'ZZJets2L2Nu muon 2012'),
        cfg('ZZJets2L2Q', dCache+'/naodell/nuTuples_v7_4/ZZJetsTo2L2Q', 5, 'ZZJets2L2Q muon 2012'),
        #cfg('ZZJets4L', dCache+'/naodell/nuTuples_v7_4/ZZJetsTo4L', 5, 'ZZJets4L muon 2012'),
        cfg('ZZTo4e', dCache+'/naodell/nuTuples_v7_4/ZZTo4e', 5, 'ZZ4e muon 2012'),
        cfg('ZZTo4mu', dCache+'/naodell/nuTuples_v7_4/ZZTo4mu', 5, 'ZZ4mu muon 2012'),
        cfg('ZZTo4tau', dCache+'/naodell/nuTuples_v7_4/ZZTo4tau', 5, 'ZZ4tau muon 2012'),
        cfg('ZZTo2e2mu', dCache+'/naodell/nuTuples_v7_4/ZZTo2e2mu', 5, 'ZZ2e2mu muon 2012'),
        cfg('ZZTo2e2tau', dCache+'/naodell/nuTuples_v7_4/ZZTo2e2tau', 5, 'ZZ2e2tau muon 2012'),
        cfg('ZZTo2mu2tau', dCache+'/naodell/nuTuples_v7_4/ZZTo2mu2tau', 5, 'ZZ2mu2tau muon 2012'),
        cfg('WWJets2L2Nu', dCache+'/naodell/nuTuples_v7_4/WWJetsTo2L2Nu', 5, 'WWJets2L2Nu muon 2012'),
        cfg('WZJets3LNu', dCache+'/naodell/nuTuples_v7_4/WZJetsTo3LNu', 5, 'WZJets3LNu muon 2012'),
        cfg('WZJets2L2Q', dCache+'/naodell/nuTuples_v7_4/WZJetsTo2L2Q', 10, 'WZJets2L2Q muon 2012'),

        cfg('QCD_20-30_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_20_30_EMEnriched', 10, 'QCD_20-30_EM muon 2012'),
        cfg('QCD_30-80_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_30_80_EMEnriched', 10, 'QCD_30-80_EM muon 2012'),
        cfg('QCD_80-170_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_80_170_EMEnriched', 10, 'QCD_80-170_EM muon 2012'),
        cfg('QCD_170-250_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_170_250_EMEnriched', 10, 'QCD_170-250_EM muon 2012'),
        cfg('QCD_250-350_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_250_350_EMEnriched', 10, 'QCD_250-350_EM muon 2012'),
        cfg('QCD_350_EM', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_350_EMEnriched', 10, 'QCD_350_EM muon 2012'),
        cfg('QCD_20_MU', dCache+'/naodell/nuTuples_v7_4/QCD_Pt_20_MuEnrichedPt_15', 10, 'QCD_20_MU muon 2012'),

        cfg('ggHToZZ4L_M-125', dCache+'/naodell/nuTuples_v7_4/GluGluToHToZZTo4L_M-125', 2, 'ggHToZZ4L_M-125 muon 2012'),
        cfg('ggHToWW2L2Nu_M-125', dCache+'/naodell/nuTuples_v7_4/GluGluToHToWWTo2LAndTau2Nu_M-125', 2, 'ggHToWW2L2Nu_M-125 muon 2012'),
        cfg('WHToWWW3L_M-125', dCache+'/naodell/nuTuples_v7_4/WH_HToWW_3l_M-125', 2, 'WHToWWW3L_M-125 muon 2012')
        ])

    signal.extend([
        cfg('FCNC_M125_tHj', dCache+'/naodell/nuTuples_v6_8TeV/FCNH_M125_t', 1, 'FCNC_M125_t mc 2012'),
        cfg('FCNC_M125_tbarHj', dCache+'/naodell/nuTuples_v6_8TeV/FCNH_M125_tbar', 1, 'FCNC_M125_tbar mc 2012')
        ])


if period == '2011':
    data.extend([
        cfg('muon_2011A', dCache+'/andreypz/nuTuples_v2_7TeV/DoubleMu_HZZ_Run2011A', 30, 'DATA muon 2011'),
        cfg('muon_2011B', dCache+'/andreypz/nuTuples_v2_7TeV/DoubleMu_HZZ_Run2011B', 30, 'DATA muon 2011'),
        cfg('muEG_2011A', dCache+'/naodell/nuTuples_v2_7TeV/MuEG_Run2011A', 30, 'DATA muEG 2011'),
        cfg('muEG_2011B', dCache+'/naodell/nuTuples_v2_7TeV/MuEG_Run2011B', 30, 'DATA muEG 2011'),
        #cfg('electron_2011A', dCache+'/naodell/nuTuples_v2_7TeV/DoubleElectron_Run2011A', 20, 'DATA 16,17,18 electron 2011')
        #cfg('electron_2011B', dCache+'/andreypz/nuTuples_v2_7TeV/DoubleElectron_HZZ_Run2011B', 20, 'DATA 16,17,18 electron 2011')
        ])

    bg.extend([
        cfg('ZJets', dCache+'/andreypz/nuTuples_v2_7TeV/DYjets', 40, 'ZJets mc 2011'),
        cfg('ttbar', dCache+'/andreypz/nuTuples_v2_7TeV/TTJets', 40, 'ttbar mc 2011'),
        cfg('tbarW', dCache+'/andreypz/nuTuples_v2_7TeV/tbarW', 30, 'tbarW mc 2011'),
        cfg('tW', dCache+'/andreypz/nuTuples_v2_7TeV/tW', 15, 'tW mc 2011'),
        cfg('ZZJets2L2Nu', dCache+'/naodell/nuTuples_v2_7TeV/ZZJetsTo2L2Nu', 5, 'ZZJets2L2Nu mc 2011'),
        cfg('ZZJets2L2Q', dCache+'/naodell/nuTuples_v2_7TeV/ZZJetsTo2L2Q', 5, 'ZZJets2L2Q mc 2011'),
        cfg('WWJets2L2Nu', dCache+'/naodell/nuTuples_v2_7TeV/WWJetsTo2L2Nu', 5, 'WWJets2L2Nu mc 2011'),
        cfg('WZJets3LNu', dCache+'/naodell/nuTuples_v2_7TeV/WZJetsTo3LNu', 5, 'WZJets3LNu mc 2011'),
        cfg('WZJets2L2Q', dCache+'/naodell/nuTuples_v2_7TeV/WZJetsTo2L2Q', 5, 'WZJets2L2Q mc 2011')
        ])


    signal.extend([
        cfg('FCNC', dCache+'/devildog/nuTuples_v2_7TeV/FCNC_tH', 5, 'FCNC mc 2011')
>>>>>>> 18ffb5e611afdfafac090643e0c149608c892737
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

