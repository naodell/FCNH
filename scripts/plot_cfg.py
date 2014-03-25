#! /usr/bin/env python
import subprocess, shlex, datetime, copy
from multiprocessing import Process

import ROOT as r
from PlotProducer import *
from TableMaker import *

now         = datetime.datetime.now()
currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

### Get command line arguements

if len(sys.argv) > 2:
    batch   = sys.argv[1]
    suffix  = sys.argv[2]
else:
    print 'Need to specify input batch and suffix!!!  Exiting...'
    exit()

### This is the config file for manipulating 
### histograms using the PlotProducer class.  

plotType    = '.png'
selection   = 'fcnh'

cutList     = ['1_preselection']
cutList.extend(['2_Z_veto', '3_2jet', '4_MET', '.'])#, '5_BDT'])
#cutList.extend(['.', '.', '.', 'X_0jet', 'X_1jet'])

crList      = []#'CR_WZ', 'CR_ttbar', 'CR_ZFake']

period      = '2012'
LUMIDATA    = 19.712 

doPlots     = True
doLog       = False
doEff       = False
doRatio     = True
doNorm      = True
do1D        = True
do2D        = True

doOS        = False
doSS        = True
do3l        = True

doYields    = True

### Categories to be plotted ###
catSS       = ['ss_inclusive']
catSS.extend(['ss_mumu', 'ss_ee', 'ss_emu'])
#catSS.extend(['ss_endcap', 'ss_mixed', 'ss_barrel'])
catOS       = ['os_inclusive']
catOS.extend(['os_mumu', 'os_ee', 'os_emu']) 
cat3l       = ['3l_inclusive']
cat3l.extend(['3l_eee', '3l_eemu', '3l_emumu', '3l_mumumu'])
#cat3l.extend(['3l_OSSF', '3l_SSSF'])

### Samples to be included in stacks ###
samples     = {'all':[], 'inclusive':[], 'os':[], 'WZ':[], 'ttbar':[], 'ttZ':[], 'ZFake':[],
                '3l_inclusive':[], '3l_eee':[], '3l_eemu':[], '3l_emumu':[], '3l_mumumu':[], 
                'ss_inclusive':[], 'ss_ee':[], 'ss_emu':[], 'ss_mumu':[]}

#samples['all'].append('higgs')
samples['all'].append('Triboson')
samples['all'].append('ttV')
samples['all'].append('Diboson')
samples['all'].append('top')
samples['all'].append('ZJets')
samples['all'].append('WG')
#samples['all'].append('ZZ4l')
#samples['all'].append('WZJets3LNu')
#samples['all'].append('QCD')
#samples['all'].extend(['ZbbToLL', 'WbbToLNu']) #, 'ZGstar'])

#samples['inclusive'].append('higgs')
samples['inclusive'].append('Triboson')
samples['inclusive'].append('ttV')
samples['inclusive'].append('Diboson')
samples['inclusive'].append('top')
samples['inclusive'].append('ZJets')

## trilepton categories
#samples['3l_inclusive'].append('higgs')
samples['3l_inclusive'].append('Triboson')
samples['3l_inclusive'].append('ttV')
samples['3l_inclusive'].append('ZZ4l')
#samples['3l_inclusive'].extend(['ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau'])
samples['3l_inclusive'].append('WZJets3LNu')

## eee
samples['3l_eee'].extend(samples['3l_inclusive'])
samples['3l_eee'].extend(['eFakes', 'llFakes'])

## eemu
samples['3l_eemu'].extend(samples['3l_inclusive'])
samples['3l_eemu'].extend(['eFakes', 'muFakes', 'llFakes'])

## emumu
samples['3l_emumu'].extend(samples['3l_inclusive'])
samples['3l_emumu'].extend(['eFakes', 'muFakes', 'llFakes'])

## mumumu
samples['3l_mumumu'].extend(samples['3l_inclusive'])
samples['3l_mumumu'].extend(['muFakes', 'llFakes'])

## inclusive
samples['3l_inclusive'].append('Fakes')

#samples['3l'].append('top')
#samples['3l'].append('ZJets')
#samples['3l'].append('Diboson')
#samples['3l'].append('ZGstar')
#samples['3l'].append('WGStar')
#samples['3l'].extend(['WGStarLNu2E', 'WGStarLNu2Mu', 'WGStarLNu2Tau'])

## same-sign categories
#samples['ss_inclusive'].append('higgs')
samples['ss_inclusive'].append('Triboson')
samples['ss_inclusive'].append('ttV')
samples['ss_inclusive'].append('ZZ4l')
#samples['ss_inclusive'].extend(['ZZ4mu', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu', 'ZZ2mu2tau', 'ZZ2e2tau'])
samples['ss_inclusive'].append('WZJets3LNu')
samples['ss_inclusive'].append('WG')

## dielectrons
samples['ss_ee'].extend(samples['ss_inclusive'])
samples['ss_ee'].append('QFlips')
samples['ss_ee'].extend(['eFakes', 'llFakes'])

## electron+muon
samples['ss_emu'].extend(samples['ss_inclusive'])
samples['ss_emu'].append('QFlips')
samples['ss_emu'].extend(['eFakes', 'muFakes', 'llFakes'])

## dimuons
samples['ss_mumu'].extend(samples['ss_inclusive'])
samples['ss_mumu'].extend(['muFakes', 'llFakes'])

## inclusive
samples['ss_inclusive'].append('Fakes')
samples['ss_inclusive'].append('QFlips')

## geometric categories
samples['ss_endcap']    = samples['ss_mumu']
samples['ss_mixed']     = samples['ss_mumu']
samples['ss_barrel']    = samples['ss_mumu']

## opposite-sign categories
samples['os'].extend(['Diboson', 'top', 'ZJets'])

samples['WZ'].extend(['WW/ZZ', 'top', 'ZJets', 'WZJets3LNu'])
samples['ttbar'].extend(['single top', 'ZJets', 'ttbar'])
samples['ttZ'].extend(['top', 'ZJets', 'ZZ4l', 'WZJets3LNu', 'ttW', 'ttG', 'ttZ'])
samples['ZFake'].extend(['ZZ4l', 'WZJets3LNu', 'Fakes'])

p_plot = []


if doPlots:

    print '\nMaking the plots...\n'

    r.gROOT.SetBatch()

    ### Initialize plot producer ###
    plotter = PlotProducer(inputFile = 'fcncAnalysis/combined_histos/{0}_cut1_{1}_{2}.root'.format(selection, period, batch), savePath = '', scale = LUMIDATA, isAFS = False)
    plotter.set_period(period)
    plotter.set_output_type(plotType)
    plotter.set_clean_fakes(True)

    ### DATASETS ###
    ### Specify the datasets you wish to stack 
    ### and overlay accordingly. 

    plotter.add_datasets(samples['all'])
    plotter._overlayList.extend(['DATA'])
    #plotter._overlayList.extend(['FCNH'])

    plotter.get_scale_factors()
    #plotter.get_scale_factors(['FCNH'])

    ### VARIABLES ###
    ### First specify the directories in which your
    ### histograms are stored.  If directories are 
    ### not used enter '' as the only entry.  Then 
    ### list all of the variable names you wish to 
    ### plot while giving a key value which is the 
    ### directory that they are located in as a key.

    plotter._directoryList1D            = ['Misc', 'Lepton', 'Lep+Jet', 'Dilepton', 'DileptonOS', 'Trilepton', 'MET', 'Jet', 'Fakes', 'GEN', '4l']
    plotter._directoryList2D            = ['2D']

    plotter._variableDict['Misc']       = ['PvMult', 'YieldByCut', 'YieldByCutRaw', 'EventWeight', 'TriggerStatus', 
                                            'FakeWeightUncertainty', 'BDT']

    plotter._variableDict['Lepton']     = ['LeptonCharge', 'LeptonFlavor', 
                                           'Lepton1Pt', 'Lepton2Pt','Lepton3Pt',
                                           'Lepton1Eta', 'Lepton2Eta', 'Lepton3Eta',
                                           'Lepton1IsoRel', 'Lepton2IsoRel', 'Lepton3IsoRel', 
                                           'ElectronPt', 'ElectronEta',
                                           'ElectronDxy', 'ElectronDz',
                                           'ElectronIsoRel', 'ElectronIso',
                                           'MuonPt', 'MuonEta',
                                           'MuonDxy', 'MuonDz', 
                                           'MuonIsoRel', 'MuonIso',
                                           'Lepton1dxy', 'Lepton1dz',
                                           'Lepton2dxy', 'Lepton2dz',
                                           'Lepton3dxy', 'Lepton3dz',
                                           'LeptonMult', 'OverlapEleMu', 'fakeableOverlapMult']
                                           #'Lepton1Phi', 'Lepton2Phi', 'Lepton3Phi']

    plotter._variableDict['Dilepton']   = ['DileptonMass21', 'DileptonTransMass21', 'DileptonQt21', 'DileptonBalance21',
                                           'DileptonDeltaPhi21', 'DileptonDeltaEta21', 'DileptonDeltaR21', 'DileptonDeltaPt21']
                                           #'DileptonMass31', 'DileptonTransMass31', 'DileptonQt31', 
                                           #'DileptonDeltaPhi31', 'DileptonDeltaEta31', 'DileptonDeltaR31', 'DileptonDeltaPt31',
                                           #'DileptonMass32', 'DileptonTransMass32', 'DileptonQt32', 
                                           #'DileptonDeltaPhi32', 'DileptonDeltaEta32', 'DileptonDeltaR32', 'DileptonDeltaPt32']

    plotter._variableDict['DileptonOS'] = ['DileptonOSMass', 'DileptonOSTransMass', 'DileptonOSBalance',
                                           'DileptonOSQt', 'DileptonOSDeltaPt', 'DileptonOSDeltaR', 
                                           'DileptonOSDeltaEta', 'DileptonOSDeltaPhi'] 

    plotter._variableDict['Trilepton']  = ['DileptonLepDeltaR', 'DileptonLepDeltaPhi', 'DileptonLepDeltaEta', 
                                           'Lep3MetMT', 'TrileptonMass', 'TrileptonPt']

    plotter._variableDict['Lep+Jet']    = ['Lepton1BJetDeltaPhi', 'Lepton1BJetDeltaEta', 'Lepton1BJetDeltaR', 'Lepton1BJetDeltaPt',
                                           #'Lepton2BJetDeltaPhi', 'Lepton2BJetDeltaEta', 'Lepton2BJetDeltaR', 'Lepton2BJetDeltaPt',
                                           #'Lepton3BJetDeltaPhi', 'Lepton3BJetDeltaEta', 'Lepton3BJetDeltaR', 'Lepton3BJetDeltaPt',
                                           'Lepton1JetDeltaPhi', 'Lepton1JetDeltaEta', 'Lepton1JetDeltaR', 'Lepton1JetDeltaPt',
                                           #'Lepton2JetDeltaPhi', 'Lepton2JetDeltaEta', 'Lepton2JetDeltaR', 'Lepton2JetDeltaPt',
                                           #'Lepton3JetDeltaPhi', 'Lepton3JetDeltaEta', 'Lepton3JetDeltaR', 'Lepton3JetDeltaPt',
                                           'DileptonJetDeltaPhi', 'DileptonJetDeltaEta', 'DileptonJetDeltaR', 'DileptonJetDeltaPt',
                                           'DileptonBJetDeltaPhi', 'DileptonBJetDeltaEta', 'DileptonBJetDeltaR', 'DileptonBJetDeltaPt',
                                           'OverlapJetMult'
                                          ]

    plotter._variableDict['Jet']        = ['Jet1Pt', 'Jet2Pt',# 'Jet3Pt',
                                           'Jet1Eta', 'Jet2Eta',# 'Jet3Eta',
                                           #'Jet1Phi', 'Jet2Phi', 'Jet3Phi',
                                           'BJet1BDiscr', 'BJet1Pt', 'BJet1Eta', #'BJet1Phi', 
                                           'BJet2BDiscr', 'BJet2Pt', 'BJet2Eta', #'BJet2Phi',
                                           'HT', 'HTs', 'EventBalance', 'Centrality',
                                           'JetBJetDeltaPhi', 'JetBJetDeltaEta', 'JetBJetDeltaR',
                                           'JetMultCharge', 'JetMult', 'BJetMult', 'AllJetMult',
                                           'MatchedMuJetBDiscr', 'MatchedEleJetBDiscr', 
                                           'DijetMass']

    plotter._variableDict['MET']        = ['Met', 'MHT', 'METLD', 'MHT-MET', 'MetPhi', 'MetSumEt',
                                           'MetLepton1DeltaPhi', 'MetLepton2DeltaPhi', 'MetLepton3DeltaPhi'
                                           'MetLepDeltaPhiMin', 'nearLepIndex', 'ProjectedMet'] 

    plotter._variableDict['Fakes']      = ['FakeablePt', 'FakeableEta', 'FakeablePhi'
                                           'FakeableDxy', 'FakeableDz', 'FakeableIsoRel']

    plotter._variableDict['GEN']        = ['GenChargeMisId', 'GenMisIdPt', 'GenMisIdEta',
                                           'GenDeltaR', 'GenBalance']

    plotter._variableDict['4l']         = ['4lMass', '4lPt', '4lSumPt', '4lMet']

    plotter._variableDict['2D']         = ['metVsHt', 'metVsSqrtHt', 'TrileptonMVsDileptonMOS',
                                            'DileptonMVsDeltaROS', 'DileptonQtVsDeltaROS',
                                            #'DileptonM13VsM21', 'DileptonM12VsM31', 'DileptonM21VsM32',
                                            #'DalitzM13VsM21', 'DalitzM12VsM31', 'DalitzM21VsM32',
                                            'LepChargeVsFlavor']


     ###################   
     ### MAKE PLOTS! ###  
     ###################   

    r.gROOT.SetStyle('Plain')
    r.gStyle.SetOptStat(0)
    #r.gROOT.ProcessLine('.L ./tdrStyle.C')
    #r.setTDRStyle()

    ### inclusive ###

    inclusive_plotter = copy.deepcopy(plotter)
    inclusive_plotter.add_datasets(samples['inclusive'], Clear=True)
    inclusive_plotter._overlayList = ['DATA'] # overlaySamples

    for i, cut in enumerate(cutList):
        if cut == '.':
            continue

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

        inclusive_plotter.make_save_path(outFile, clean=True)
        p_plot.append(Process(name = cut[2:] + '/inclusive', target = plotter_wrapper, args=(inclusive_plotter, 'inclusive', inFile, outFile, do1D, do2D, False, doLog, doRatio, False)))


    ### 3l selection ###
    if do3l:
        for category in cat3l:
            plotter_3l = copy.deepcopy(plotter)
            plotter_3l.add_datasets(samples[category], Clear=True)
            plotter_3l._overlayList = ['DATA']#, 'FCNH']

            for i, cut in enumerate(cutList):
                if cut == '.':
                    continue

                inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

                if doLog:
                    outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
                else:
                    outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

                plotter_3l.make_save_path(outFile, clean=True)
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(plotter_3l, category, inFile, outFile, do1D, do2D, False, doLog, doRatio, doEff)))

    ### ss selection ###
    if doSS:
        for category in catSS:
            ss_plotter = copy.deepcopy(plotter)
            ss_plotter.add_datasets(samples[category], Clear=True)
            ss_plotter._overlayList = ['DATA']#, 'FCNH']
            #ss_plotter.set_clean_fakes(False)

            for i, cut in enumerate(cutList):
                if cut == '.':
                    continue

                inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, i+1, period, batch)

                if doLog:
                    outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
                else:
                    outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

                ss_plotter.make_save_path(outFile, clean=True)
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(ss_plotter, category, inFile, outFile, do1D, do2D, False, doLog, doRatio, doEff)))

                ### overlay of fake distributions for single and double fakes ###
                if cut in ['1_preselection', 'X_0jet', 'X_1jet', '3_2jet']:
                    fake_plotter = copy.deepcopy(plotter)
                    fake_plotter._overlayList = ['llFakes', 'muFakes']
                    fake_plotter.draw_normalized(True)

                    inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, i+1, period, batch)
                    if doLog:
                        outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix,  'fakes_overlay_{0}'.format(cut[2:]))
                    else:
                        outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'fakes_overlay_{0}'.format(cut[2:]))
                    fake_plotter.make_save_path(outFile, clean=True)

                    p_plot.append(Process(name = 'fakes_overlay_{0}/{1}'.format(cut[2:], category), target = plotter_wrapper, args=(fake_plotter, category, inFile, outFile, do1D, False, True, True, False, False)))


    ### os selection ###
    if doOS:
        os_plotter = copy.deepcopy(plotter)
        os_plotter.add_datasets(samples['os'], Clear=True)
        os_plotter._overlayList = ['DATA']

        for i, cut in enumerate(cutList):
            if cut == '.':
                continue

            inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

            if doLog:
                outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
            else:
                outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

            os_plotter.make_save_path(outFile, clean=True)

            for category in catOS:
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(os_plotter, category, inFile, outFile, do1D, do2D, False, doLog, doRatio, False)))

    doLog = False

    ### WZ control region
    if 'CR_WZ' in crList:
        wz_plotter = copy.deepcopy(plotter)
        wz_plotter.add_datasets(samples['WZ'], Clear=True)
        wz_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 7, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'CR_WZ')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'CR_WZ')

        wz_plotter.make_save_path(outFile, clean=True)

        for category in cat3l:
            p_plot.append(Process(name = 'CR_WZ/' + category, target = plotter_wrapper, args=(wz_plotter, category, inFile, outFile, do1D, False, False, doLog, doRatio, False)))

    ### ttbar control region
    if 'CR_ttbar' in crList:
        ttbar_plotter = copy.deepcopy(plotter)
        ttbar_plotter.add_datasets(samples['ttbar'],  Clear=True)
        ttbar_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 8, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'CR_ttbar')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'CR_ttbar')

        ttbar_plotter.make_save_path(outFile, clean=True)

        p_plot.append(Process(name = 'CR_ttbar/os_emu', target = plotter_wrapper, args=(ttbar_plotter, 'os_emu', inFile, outFile, do1D, False, False, doLog, doRatio, False)))

    doLog = True

    ### ZPlusFake control region
    if 'CR_ZFake' in crList:
        ZFake_plotter = copy.deepcopy(plotter)
        ZFake_plotter.add_datasets(samples['ZFake'],  Clear=True)
        ZFake_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 9, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'CR_ZFake')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'CR_ZFake')

        ZFake_plotter.make_save_path(outFile, clean=True)

        for category in cat3l:
            p_plot.append(Process(name = 'CR_ZFake/' + category, target = plotter_wrapper, args=(ZFake_plotter, category, inFile, outFile, do1D, False, False, doLog, doRatio, False)))
    
    ### low delta eta ss control region
    if 'high_mass_ss' in crList:
        hm_plotter = copy.deepcopy(plotter)
        hm_plotter.add_datasets(samples['ss'],  Clear=True)
        hm_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 9, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'high_mass_ss')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'high_mass_ss')

        hm_plotter.make_save_path(outFile, clean=True)

        for category in catSS:
            p_plot.append(Process(name = 'high_mass_ss/' + category, target = plotter_wrapper, args=(hm_plotter, category, inFile, outFile, do1D, False, doLog, doRatio, False)))

    ### low delta eta ss control region
    if 'low_mass_ss' in crList:
        lm_plotter = copy.deepcopy(plotter)
        lm_plotter.add_datasets(samples['ss'],  Clear=True)
        lm_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 10, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'low_mass_ss')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'low_mass_ss')

        lm_plotter.make_save_path(outFile, clean=True)

        for category in catSS:
            p_plot.append(Process(name = 'high_mass_ss/' + category, target = plotter_wrapper, args=(lm_plotter, category, inFile, outFile, do1D, False, doLog, doRatio, False)))

    ### low delta eta ss control region
    if 'barrel_leptons' in crList:
        bl_plotter = copy.deepcopy(plotter)
        bl_plotter.add_datasets(samples['ss'],  Clear=True)
        bl_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 11, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'barrel_leptons')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'barrel_leptons')

        bl_plotter.make_save_path(outFile, clean=True)


### End of configuration for PlotProducer ###

for process in p_plot:
    print 'Plotting {0}'.format(process.name)
    process.start()

for process in p_plot:
    process.join()

print '\n'

if doYields:
    ### Initialize table maker ###
    tableFile       = file('yields/.yields_tmp.tex', 'w')
    yieldTable      = TableMaker('fcncAnalysis/combined_histos/{0}_cut1_{1}_{2}.root'.format(selection, period, batch), tableFile, scale = LUMIDATA, delimiter = '&', doSumBG = True)
    yieldTable.set_period(period)

    yieldTable.add_datasets(samples['all'], Clear = True)

    if not doPlots:
        yieldTable.get_scale_factors()
        #yieldTable.get_scale_factors(['FCNH'])

    if do3l:
        #yieldTable._columnList  = ['Irreducible', 'Fakes', 'BG', 'DATA', 'FCNH']#, 'Significance'] 
        yieldTable._columnList  = samples['3l_inclusive'] + ['BG', 'DATA']#, 'FCNH']#, 'Significance'] 
        #yieldTable._columnList  = ['BG', 'DATA', 'FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_M125_t_semilep', 'FCNC_M125_t_ZZ', 'FCNC_M125_t_TauTau','FCNH']# 'Significance'] 
        #yieldTable._columnList  = ['BG', 'DATA', 'FCNH']#, 'Significance'] 

        #yieldTable.add_datasets(['Irreducible', 'Fakes'], Clear = True)
        yieldTable.add_datasets(samples['3l_inclusive'], Clear = True)
        #yieldTable.add_datasets('FCNH')
        yieldTable.add_datasets('DATA')

        yieldTable._rowList = 5*['.'] + ['ss dilepton', 'Z removal', '2+ jets'] + 6*['.'] + ['0-jet', '1-jet']# + 7*['.'] + ['BDT']

        for category in cat3l:
            yieldTable._category = category
            histDict = yieldTable.get_hist_dict('YieldByCut')
            yieldTable.print_table(histDict, doErrors = True, doEff = False, startBin = 1)

    if doSS:
        #yieldTable._columnList  = ['Irreducible', 'Fakes', 'QFlips', 'BG', 'DATA', 'FCNH']#, 'Significance'] 
        yieldTable._columnList  = samples['ss_inclusive'] + ['BG', 'DATA']#, 'FCNH']#, 'Significance'] 
        #yieldTable._columnList  = ['BG', 'DATA', 'FCNH']#, 'Significance'] 

        #yieldTable.add_datasets(['Irreducible', 'Fakes', 'QFlips'], Clear = True)
        yieldTable.add_datasets(samples['ss_inclusive'], Clear = True)
        #yieldTable.add_datasets('FCNH')
        yieldTable.add_datasets('DATA')
        yieldTable._rowList = 5*['.'] + ['ss dilepton', 'Z removal', '2+ jets', 'MET'] + 5*['.'] + ['0-jet', '1-jet']# + 7*['.'] + ['BDT']

        for category in catSS:
            yieldTable._category = category
            histDict = yieldTable.get_hist_dict('YieldByCut')
            yieldTable.print_table(histDict, doErrors = True, doEff = False, startBin = 1)

    crCats = {'CR_WZ':'3l_inclusive', 'CR_ttbar':'os_emu', 'CR_ttZ':'3l_inclusive'}
    for i,CR in enumerate(crList):
        if CR[:2] != 'CR' or CR == 'CR_ZFake': continue

        yieldTable.set_input_file('fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, i+5, period, batch))
        yieldTable._columnList  = samples[CR[3:]] + ['BG', 'DATA']

        yieldTable.add_datasets(samples[CR[3:]], Clear = True)
        yieldTable.add_datasets('DATA')

        yieldTable._rowList = ['preselection'] + (4+i)*['.'] + [CR[3:]]

        yieldTable._category = crCats[CR]
        histDict = yieldTable.get_hist_dict('YieldByCut')
        yieldTable.print_table(histDict, doErrors = True, doEff = False, startBin = 6)

    ### Special case for ZZ->4l control region ###
    #yieldTable.set_input_file('fcncAnalysis/combined_histos/{0}_cut1_{1}_{2}.root'.format(selection, period, batch))
    #yieldTable.add_datasets(['ZZ4l', 'DATA'], Clear = True)
    #yieldTable._columnList  = ['ZZ4l'] + ['BG', 'DATA']
    #yieldTable._rowList = 8*['.'] + ['ZZ4l']

    #yieldTable._category = 'inclusive'
    #histDict = yieldTable.get_hist_dict('YieldByCut')
    #yieldTable.print_table(histDict, doErrors = False, doEff = False, startBin = 6)

    tableFile.close()

    subprocess.call('pdflatex -output-dir=yields yields/yields.tex', shell = True)
    subprocess.call('cp yields/yields.pdf plots/{0}/{1}_{2}_{3}/.'.format(currentDate, selection, batch, suffix), shell = True)
    subprocess.call('cp yields/.yields_tmp.tex plots/{0}/{1}_{2}_{3}/yields.tex'.format(currentDate, selection, batch, suffix), shell = True)
