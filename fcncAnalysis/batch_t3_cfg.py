#! /usr/bin/env python
import BatchMaster as b
import sys

cfg = b.JobConfig

''' Specify parameters '''
dataDir     = '/tthome/naodell/storage/ntuples/Data'
signalDir   = '/tthome/naodell/storage/ntuples/MC_FCNH/'
mcDir       = '/tthome/naodell/storage/data/nuTuples_v9.6_8TeV/MC/'
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
        cfg('muon_2012A', '{0}/DoubleMu_Run2012A_v2'.format(dataDir), 40, 'DATA_MUON muon 2012'),
        cfg('muon_2012B', '{0}/DoubleMu_Run2012B'.format(dataDir), 40, 'DATA_MUON muon 2012'),
        cfg('muon_2012C', '{0}/DoubleMu_Run2012C'.format(dataDir), 40, 'DATA_MUON muon 2012'),
        cfg('muon_2012D', '{0}/DoubleMu_Run2012D'.format(dataDir), 45, 'DATA_MUON muon 2012'),

        cfg('electron_2012A', '{0}/DoubleElectron_Run2012A'.format(dataDir), 40, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012B', '{0}/DoubleElectron_Run2012B'.format(dataDir), 40, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012C', '{0}/DoubleElectron_Run2012C'.format(dataDir), 40, 'DATA_ELECTRON electron 2012'),
        cfg('electron_2012D', '{0}/DoubleElectron_Run2012D'.format(dataDir), 45, 'DATA_ELECTRON electron 2012'),

        cfg('muEG_2012A', '{0}/MuEG_Run2012A'.format(dataDir), 40, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012B', '{0}/MuEG_Run2012B'.format(dataDir), 40, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012C', '{0}/MuEG_Run2012C'.format(dataDir), 40, 'DATA_MUEG muEG 2012'),
        cfg('muEG_2012D', '{0}/MuEG_Run2012D'.format(dataDir), 45, 'DATA_MUEG muEG 2012')
        ])

    bg.extend([
        #cfg('ZJetsToTauTau_M-50',      '{0}/DYJetsToLL_M-50'.format(mcDir),            60,  'DYJetsToTauTau      {0}   2012'.format(mcTrigger)),
        cfg('ZJets_M-50',      '{0}/DYJetsToLL_M-50'.format(mcDir),            30,  'ZJets_M-50      {0}   2012'.format(mcTrigger)),
        cfg('ZJets_M-10To50',  '{0}/DYJetsToLL_M-10To50filter'.format(mcDir),  20,  'ZJets_M-10To50  {0}   2012'.format(mcTrigger)),
        #cfg('WbbToLNu',       '{0}/WbbJetsToLNu'.format(mcDir),               20,  'WbbToLNu        {0}   2012'.format(mcTrigger)),
        #cfg('WjetToLNu',       '{0}/WJetsToLNu'.format(mcDir),                 50,  'WJetsToLNu      {0}   2012'.format(mcTrigger)),
        #cfg('W1JetToLNu',        '{0}/W1JetsToLNu'.format(mcDir),                 10,  'W1JetsToLNu      {0}   2012'.format(mcTrigger)),
        #cfg('W2JetToLNu',        '{0}/W2JetsToLNu'.format(mcDir),                 10,  'W2JetsToLNu      {0}   2012'.format(mcTrigger)),
        #cfg('W3JetToLNu',        '{0}/W3JetsToLNu'.format(mcDir),                 10,  'W3JetsToLNu      {0}   2012'.format(mcTrigger)),
        #cfg('W4JetToLNu',        '{0}/W4JetsToLNu'.format(mcDir),                 10,  'W4JetsToLNu      {0}   2012'.format(mcTrigger)),
        #cfg('WGStarLNu2Mu',    '{0}/WGstarToLNu2Mu'.format(mcDir),             5,   'WGStarLNu2Mu    {0}   2012'.format(mcTrigger)),
        #cfg('WGStarLNu2Tau',   '{0}/WGstarToLNu2Tau'.format(mcDir),            5,   'WGStarLNu2Tau   {0}   2012'.format(mcTrigger)),

        #cfg('WGStarLNu2E', mcDir+'/WGstarToLNu2E', 5, 'WGStarLNu2E muon 2012'),
        #cfg('ZbbToLL', mcDir+'/ZbbToLL', 20, 'ZbbToLL muon 2012'),
        #cfg('ZG', mcDir+'/ZGToLLG', 10, 'ZG muon 2012'),
        #cfg('WG', '{0}/WGToLNuG'.format('/tthome/naodell/storage/ntuples/MC_skimmed'), 10, 'WG {0} 2012'.format(mcTrigger)),

        #cfg('ttbarHad',        '{0}/TTJets'.format(mcDir),   30,  'ttbarHad        {0}  2012'.format(mcTrigger)),
        #cfg('ttbarLep',        '{0}/TTJets'.format(mcDir),   30,  'ttbarLep        {0}  2012'.format(mcTrigger)),
        #cfg('tbarW',           '{0}/Tbar_tW'.format(mcDir),  10,  'tbarW           {0}  2012'.format(mcTrigger)),
        #cfg('tW',              '{0}/T_tW'.format(mcDir),     10,  'tW              {0}  2012'.format(mcTrigger)),
        cfg('ttW',             '{0}/TTWJets'.format(mcDir),  5,   'ttW             {0}  2012'.format(mcTrigger)),
        cfg('ttZ',             '{0}/TTZJets'.format(mcDir),  5,   'ttZ             {0}  2012'.format(mcTrigger)),
        cfg('ttG',             '{0}/TTGJets'.format(mcDir),  5,   'ttG             {0}  2012'.format(mcTrigger)),
        #cfg('t_t-channel',     '{0}/T_t'.format(mcDir),      5,   't_t-channel     {0}  2012'.format(mcTrigger)),
        #cfg('tbar_t-channel',  '{0}/Tbar_t'.format(mcDir),   5,   'tbar_t-channel  {0}  2012'.format(mcTrigger)),

        cfg('WWW',  '{0}/WWWJets'.format(mcDir),         5,  'WWW  {0}  2012'.format(mcTrigger)),
        cfg('WWZ',  '{0}/WWZNoGstarJets'.format(mcDir),  5,  'WWZ  {0}  2012'.format(mcTrigger)),
        cfg('WZZ',  '{0}/WZZNoGstarJets'.format(mcDir),  5,  'WZZ  {0}  2012'.format(mcTrigger)),
        cfg('ZZZ',  '{0}/ZZZNoGstarJets'.format(mcDir),  5,  'ZZZ  {0}  2012'.format(mcTrigger)),
        cfg('WWG',  '{0}/WWGJets'.format(mcDir),         5,  'WWG  {0}  2012'.format(mcTrigger)),

        #cfg('WWJets2L2Nu',  '{0}/WWJetsTo2L2Nu'.format(mcDir),  10,  'WWJets2L2Nu  {0}  2012'.format(mcTrigger)),
        #cfg('ZZJets2L2Nu',  '{0}/ZZJetsTo2L2Nu'.format(mcDir),  10,  'ZZJets2L2Nu  {0}  2012'.format(mcTrigger)),
        #cfg('WZJets2L2Q',   '{0}/WZJetsTo2L2Q'.format(mcDir),   10,  'WZJets2L2Q   {0}  2012'.format(mcTrigger)),
        #cfg('ZZJets2L2Q',   '{0}/ZZJetsTo2L2Q'.format(mcDir),   10,  'ZZJets2L2Q   {0}  2012'.format(mcTrigger)),
        cfg('WZJets3LNu',   '{0}/WZJetsTo3LNu'.format(mcDir),   10,  'WZJets3LNu   {0}  2012'.format(mcTrigger)),
        cfg('ZZTo4e',       '{0}/ZZTo4e'.format(mcDir),         10,  'ZZ4e         {0}  2012'.format(mcTrigger)),
        cfg('ZZTo4mu',      '{0}/ZZTo4mu'.format(mcDir),        10,  'ZZ4mu        {0}  2012'.format(mcTrigger)),
        cfg('ZZTo4tau',     '{0}/ZZTo4tau'.format(mcDir),       10,  'ZZ4tau       {0}  2012'.format(mcTrigger)),
        cfg('ZZTo2e2mu',    '{0}/ZZTo2e2mu'.format(mcDir),      10,  'ZZ2e2mu      {0}  2012'.format(mcTrigger)),
        cfg('ZZTo2e2tau',   '{0}/ZZTo2e2tau'.format(mcDir),     10,  'ZZ2e2tau     {0}  2012'.format(mcTrigger)),
        cfg('ZZTo2mu2tau',  '{0}/ZZTo2mu2tau'.format(mcDir),    10,  'ZZ2mu2tau    {0}  2012'.format(mcTrigger)),

        cfg('WmWmqq',  '{0}/WmWmqq'.format(mcDir),    10,  'WmWmqq    {0}  2012'.format(mcTrigger)),
        cfg('WpWpqq',  '{0}/WpWpqq'.format(mcDir),    10,  'WpWpqq    {0}  2012'.format(mcTrigger)),
        cfg('WWDPS',   '{0}/WWDPS'.format(mcDir),     10,  'WWDPS  {0}  2012'.format(mcTrigger)),
        cfg('TBZ',     '{0}/TBZ'.format(mcDir),       10,  'TBZ {0} 2012'.format(mcTrigger)),

        #cfg('QCD_20-30_EM',    '{0}/QCD_Pt_20_30_EMEnriched'.format(mcDir),    40,  'QCD_20-30_EM    {0}  2012'.format(mcTrigger)),
        #cfg('QCD_30-80_EM',    '{0}/QCD_Pt_30_80_EMEnriched'.format(mcDir),    40,  'QCD_30-80_EM    {0}  2012'.format(mcTrigger)),
        #cfg('QCD_80-170_EM',   '{0}/QCD_Pt_80_170_EMEnriched'.format(mcDir),   40,  'QCD_80-170_EM   {0}  2012'.format(mcTrigger)),
        #cfg('QCD_170-250_EM',  '{0}/QCD_Pt_170_250_EMEnriched'.format(mcDir),  40,  'QCD_170-250_EM  {0}  2012'.format(mcTrigger)),
        #cfg('QCD_250-350_EM',  '{0}/QCD_Pt_250_350_EMEnriched'.format(mcDir),  40,  'QCD_250-350_EM  {0}  2012'.format(mcTrigger)),
        #cfg('QCD_350_EM',      '{0}/QCD_Pt_350_EMEnriched'.format(mcDir),      40,  'QCD_350_EM      {0}  2012'.format(mcTrigger)),
        #cfg('QCD_20_MU',       '{0}/QCD_Pt_20_MuEnrichedPt_15'.format(mcDir),  40,  'QCD_20_MU       {0}  2012'.format(mcTrigger)),

        #cfg('ggHToZZ4L_M-125', mcDir+'/GluGluToHToZZTo4L_M-125', 5, 'ggHToZZ4L_M-125 muon 2012'),
        #cfg('ggHToWW2L2Nu_M-125', mcDir+'/GluGluToHToWWTo2LAndTau2Nu_M-125', 5, 'ggHToWW2L2Nu_M-125 muon 2012'),
        #cfg('WHToWWW3L_M-125', mcDir+'/WH_HToWW_3l_M-125', 5, 'WHToWWW3L_M-125 muon 2012'),
        #cfg('TTH_M-125', mcDir+'/TTH_Inclusive_M-125', 5, 'TTH_M-125 muon 2012')

        ])

    signal.extend([
        cfg('FCNC_M125_tHj', '{0}/FCNHToWW_semileptonic_Plus'.format(signalDir), 1, 'FCNC_M125_t mc 2012'),
        cfg('FCNC_M125_tbarHj', '{0}/FCNHToWW_semileptonic_Minus'.format(signalDir), 1, 'FCNC_M125_tbar mc 2012'),
        cfg('FCNC_M125_ZZ_Plus', '{0}/FCNHToZZ_Plus'.format(signalDir), 1, 'FCNC_ZZ_t mc 2012'),
        cfg('FCNC_M125_ZZ_Minus', '{0}/FCNHToZZ_Minus'.format(signalDir), 1, 'FCNC_ZZ_tbar mc 2012'),
        cfg('FCNC_M125_TauTau_Plus', '{0}/FCNHToTauTau_Plus'.format(signalDir), 1, 'FCNC_TauTau_t mc 2012'),
        cfg('FCNC_M125_TauTau_Minus', '{0}/FCNHToTauTau_Minus'.format(signalDir), 1, 'FCNC_TauTau_tbar mc 2012'),
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

