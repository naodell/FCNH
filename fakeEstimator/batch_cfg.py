#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
dCache      = '/pnfs/cms/WAX/11/store/user'
outputPath  = '/eos/uscms/store/user/naodell/BatchOutput/fakes'
executable  = 'execBatch.csh'

selection = 'fakes'
period    = '2012'
doData    = False
doBG      = False
doSignal  = False

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

data = []
bg = []
signal = []

if period == '2012':
    data.extend([
        cfg('muon_2012A', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012A_v2', 20, 'DATA muon 2012'),
        cfg('muon_2012B', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012B_v2', 20, 'DATA muon 2012'),
        cfg('muon_2012C_prompt', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_Prompt_v2', 20, 'DATA muon 2012'),
        cfg('muon_2012C_recover', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_Recover_v2', 5, 'DATA muon 2012'),
        cfg('muon_2012C_24Aug', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012C_24Aug_v2', 10, 'DATA muon 2012'),
        cfg('muon_2012D', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleMu_Run2012D_v2', 20, 'DATA muon 2012'),

        cfg('electron_2012A', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012A_v2', 20, 'DATA electron 2012'),
        cfg('electron_2012B', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012B_v2', 20, 'DATA electron 2012'),
        cfg('electron_2012C_prompt', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_Prompt_v2', 20, 'DATA electron 2012'),
        cfg('electron_2012C_recover', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_Recover_v2', 5, 'DATA electron 2012'),
        cfg('electron_2012C_24Aug', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012C_24Aug_v2', 10, 'DATA electron 2012')
        #cfg('electron_2012D', dCache+'/naodell/nuTuples_v5_5_8TeV/DoubleElectron_Run2012D_v2', 20, 'DATA electron 2012'),

        #cfg('muEG_2012A', dCache+'/naodell/nuTuples_v5_5_8TeV/MuEG_Run2012A', 20, 'DATA muEG 2012'),
        #cfg('muEG_2012B', dCache+'/naodell/nuTuples_v5_5_8TeV/MuEG_Run2012B', 20, 'DATA muEG 2012'),
        #cfg('muEG_2012C_prompt', dCache+'/naodell/nuTuples_v5_5_8TeV/MuEG_Run2012C_prompt', 20, 'DATA muEG 2012'),
        #cfg('muEG_2012C_recover', dCache+'/naodell/nuTuples_v5_5_8TeV/MuEG_Run2012C_recover', 5, 'DATA muEG 2012'),
        #cfg('muEG_2012C_24Aug', dCache+'/naodell/nuTuples_v5_5_8TeV/MuEG_Run2012C_24Aug', 10, 'DATA muEG 2012'),
        #cfg('muEG_2012D', dCache+'/naodell/nuTuples_v5_5_8TeV/MuEG_Run2012D', 20, 'DATA muEG 2012')
        ])

    bg.extend([
        cfg('ZJets', dCache+'/andreypz/nuTuples_v5_8TeV/DYjets', 40, 'ZJets muon 2012'),
        cfg('ZJets_M-10To50', dCache+'/naodell/nuTuples_v5_8TeV/DYJetsToLL_M-10To50', 20, 'ZJets_M-10To50 muon 2012'),
        #cfg('WJets', dCache+'/naodell/nuTuples_v5_8TeV/WJetsToLNu', 30, 'WJets muon 2012'),
        #cfg('WG', dCache+'/naodell/nuTuples_v5_5_8TeV/WGToLNuG', 10, 'WG muon 2012'),
        #cfg('ZG', dCache+'/naodell/nuTuples_v5_5_8TeV/ZGToLLG', 10, 'ZG muon 2012'),
        cfg('ttbar', dCache+'/andreypz/nuTuples_v5_8TeV/TTJets', 40, 'ttbar muon 2012'),
        cfg('ttW', dCache+'/naodell/nuTuples_v5_5_8TeV/TTWJets', 10, 'ttW muon 2012'),
        cfg('ttZ', dCache+'/naodell/nuTuples_v5_5_8TeV/TTZJets', 10, 'ttZ muon 2012'),
        cfg('tbarW', dCache+'/andreypz/nuTuples_v5_8TeV/tbarW', 30, 'tbarW muon 2012'),
        cfg('tW', dCache+'/andreypz/nuTuples_v5_8TeV/tW', 15, 'tW muon 2012'),
        cfg('ZZJets2L2Nu', dCache+'/andreypz/nuTuples_v5_8TeV/ZZJetsTo2L2Nu', 5, 'ZZJets2L2Nu muon 2012'),
        cfg('ZZJets2L2Q', dCache+'/naodell/nuTuples_v5_8TeV/ZZJetsTo2L2Q', 5, 'ZZJets2L2Q muon 2012'),
        cfg('ZZJets4L', dCache+'/naodell/nuTuples_v5_8TeV/ZZJetsTo4L', 5, 'ZZJets4L muon 2012'),
        cfg('WWJets2L2Nu', dCache+'/andreypz/nuTuples_v5_8TeV/WWJetsTo2L2Nu', 5, 'WWJets2L2Nu muon 2012'),
        cfg('WZJets3LNu', dCache+'/andreypz/nuTuples_v5_8TeV/WZJetsTo3LNu', 5, 'WZJets3LNu muon 2012'),
        cfg('WZJets2L2Q', dCache+'/naodell/nuTuples_v5_8TeV/WZJetsTo2L2Q', 10, 'WZJets2L2Q muon 2012')
        #cfg('QCD_EM', dCache+'/naodell/nuTuples_v5_5_8TeV/QCD_Pt-40_doubleEMEnriched', 10, 'QCD_EM mc 2012')
        ])


    signal.extend([
        #cfg('FCNC_M125', dCache+'/devildog/nuTuples_v2_7TeV/FCNC_tH', 5, 'FCNC_M125 mc 2012')
        cfg('FCNC_M125', dCache+'/naodell/nuTuples_v5_5_8TeV/FCNC_M125_8TeV', 5, 'FCNC_M125 mc 2012'),
        cfg('FCNC_M145', dCache+'/naodell/nuTuples_v5_5_8TeV/FCNC_M145_8TeV', 5, 'FCNC_M145 mc 2012')
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
        ])

inputSamples = []

if doData:
    inputSamples.extend(data)
if doBG:
    inputSamples.extend(bg)
if doSignal:
    inputSamples.extend(signal)

if len(inputSamples) is not 0:
    batcher = b.BatchMaster(inputSamples, outputPath, shortQueue = False, stageDir = 'batchStage/outputBatch', executable = executable, selection = selection + '_' + period)
    batcher.submit_to_batch()

