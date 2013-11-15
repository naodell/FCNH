#! /usr/bin/env python
import subprocess, shlex, datetime, copy
from multiprocessing import Process

import ROOT as r
from PlotProducer import *
from TableMaker import *

now         = datetime.datetime.now()
currentDate = '{0:02d}/{1:02d}/{2:02d}'.format(now.year, now.month, now.day)

### Get command line arguements

if len(sys.argv) > 1:
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
cutList.extend(['2_Z_veto', '3_MET', '4_HT', '5_bjet'])
crList      = ['CR_WZ', 'CR_ttbar', 'CR_ttZ']

period      = '2012'
LUMIDATA    = 19.712 

doPlots     = True
doLog       = True
doEff       = False
doRatio     = True
do1D        = True
do2D        = True

doOS        = True
doSS        = True
do3l        = True

doYields    = True

### Categories to be plotted ###
catSS       = ['ss_inclusive']
catSS.extend(['ss_mumu', 'ss_ee', 'ss_emu'])
catOS       = ['os_inclusive']
catOS.extend(['os_mumu', 'os_ee', 'os_emu']) 
cat3l       = ['3l_inclusive']
cat3l.extend(['3l_OSSF', '3l_SSSF'])
cat3l.extend(['3l_eee', '3l_eemu', '3l_emumu', '3l_mumumu'])

### Samples to be included in stacks ###
samples     = {'all':[], 'inclusive':[], '3l':[], 'ss':[], 'os':[], 'WZ':[], 'ttbar':[], 'ttZ':[]}

samples['all'].append('higgs')
samples['all'].append('Triboson')
samples['all'].append('ttV')
samples['all'].append('Diboson')
samples['all'].append('top')
samples['all'].append('ZJets')
samples['all'].append('QCD')
#samples['all'].extend(['ZbbToLL', 'WbbToLNu']) #, 'ZGstar'])

samples['inclusive'].append('higgs')
samples['inclusive'].append('Triboson')
samples['inclusive'].append('ttV')
samples['inclusive'].append('Diboson')
samples['inclusive'].append('top')
samples['inclusive'].append('ZJets')

samples['3l'].append('higgs')
samples['3l'].append('Triboson')
samples['3l'].append('ttV')
samples['3l'].append('ZZ4l')
samples['3l'].append('WZJets3LNu')
samples['3l'].append('Fakes')

#samples['3l'].append('top')
#samples['3l'].append('ZJets')
#samples['3l'].append('Diboson')
#samples['3l'].append('ZGstar')
#samples['3l'].append('WGStar')
#samples['3l'].extend(['WGStarLNu2E', 'WGStarLNu2Mu', 'WGStarLNu2Tau'])

samples['ss'].append('higgs')
samples['ss'].append('ttV')
samples['ss'].append('ZZ4l')
samples['ss'].append('WZJets3LNu')
samples['ss'].append('Fakes')
samples['ss'].append('QFlips')

#samples['ss'].append('Diboson')
#samples['ss'].append('top')
#samples['ss'].append('ZJets')
#samples['ss'].append('QCD')
#samples['ss'].append('QCD_EM')
#samples['ss'].append('QCD_20_MU')
#samples['ss'].extend(['ZbbToLL', 'WbbToLNu'])

samples['os'].extend(['Diboson', 'top', 'ZJets'])

samples['WZ'].extend(['WW/ZZ', 'top', 'ZJets', 'WZJets3LNu'])
samples['ttbar'].extend(['single top', 'ZJets', 'ttbar'])
samples['ttZ'].extend(['top', 'ZJets', 'WZJets3LNu', 'ttW', 'ttG', 'ttZ'])

p_plot = []


if doPlots:

    print '\nMaking the plots...\n'

    r.gROOT.SetBatch()

    ### Initialize plot producer ###
    plotter = PlotProducer(inputFile = 'fcncAnalysis/combined_histos/{0}_cut1_{1}_{2}.root'.format(selection, period, batch), savePath = '', scale = LUMIDATA, isAFS = False)
    plotter.set_period(period)
    plotter.set_output_type(plotType)

    ### DATASETS ###
    ### Specify the datasets you wish to stack 
    ### and overlay accordingly. 

    plotter.add_datasets(samples['all'])
    plotter._overlayList.extend(['DATA'])
    plotter._overlayList.extend(['FCNH'])

    plotter.get_scale_factors(['FCNH'])

    #plotter.get_scale_factors()

    ### VARIABLES ###
    ### First specify the directories in which your
    ### histograms are stored.  If directories are 
    ### not used enter '' as the only entry.  Then 
    ### list all of the variable names you wish to 
    ### plot while giving a key value which is the 
    ### directory that they are located in as a key.

    plotter._directoryList1D            = ['Misc', 'Lepton', 'Lep+Jet', 'Dilepton', 'DileptonOS', 'Trilepton', 'MET', 'Jet', 'GEN', '4l']
    plotter._directoryList2D            = ['2D']

    plotter._variableDict['Misc']       = ['PvMult', 'YieldByCut', 'YieldByCutRaw', 'EventWeight', 'TriggerStatus', 
                                            'BDT']

    plotter._variableDict['Lepton']     = ['LeptonCharge', 'LeptonFlavor', 
                                           'Lepton1Pt', 'Lepton2Pt','Lepton3Pt',
                                           'Lepton1Eta', 'Lepton2Eta', 'Lepton3Eta',
                                           'ElectronPt', 'ElectronEta',
                                           'MuonPt', 'MuonEta',
                                           'Lepton1 dxy', 'Lepton1 dz',
                                           'Lepton2 dxy', 'Lepton2 dz',
                                           'Lepton3 dxy', 'Lepton3 dz',
                                           'LeptonMult']
                                           #'Lepton1Phi', 'Lepton2Phi', 'Lepton3Phi']

    plotter._variableDict['Dilepton']   = ['DileptonMass21', 'DileptonTransMass21', 'DileptonQt21',
                                           'DileptonDeltaPhi21', 'DileptonDeltaEta21', 'DileptonDeltaR21', 'DileptonDeltaPt21',
                                           'DileptonMass31', 'DileptonTransMass31', 'DileptonQt31', 
                                           'DileptonDeltaPhi31', 'DileptonDeltaEta31', 'DileptonDeltaR31', 'DileptonDeltaPt31',
                                           'DileptonMass32', 'DileptonTransMass32', 'DileptonQt32', 
                                           'DileptonDeltaPhi32', 'DileptonDeltaEta32', 'DileptonDeltaR32', 'DileptonDeltaPt32']

    plotter._variableDict['DileptonOS'] = ['DileptonOSMass', 'DileptonOSTransMass', 'DileptonOSBalance',
                                           'DileptonOSQt', 'DileptonOSDeltaPt', 'DileptonOSDeltaR', 
                                           'DileptonOSDeltaEta', 'DileptonOSDeltaPhi'] 

    plotter._variableDict['Trilepton']  = ['DileptonLepDeltaR', 'DileptonLepDeltaPhi', 'DileptonLepDeltaEta', 'Lep3MetMT', 'TrileptonMass']

    plotter._variableDict['Lep+Jet']    = ['Lepton1BJetDeltaPhi', 'Lepton1BJetDeltaEta', 'Lepton1BJetDeltaR', 'Lepton1BJetDeltaPt',
                                           'Lepton2BJetDeltaPhi', 'Lepton2BJetDeltaEta', 'Lepton2BJetDeltaR', 'Lepton2BJetDeltaPt',
                                           'Lepton3BJetDeltaPhi', 'Lepton3BJetDeltaEta', 'Lepton3BJetDeltaR', 'Lepton3BJetDeltaPt'
                                           ]


    plotter._variableDict['Jet']        = ['Jet1Pt', 'Jet2Pt',# 'Jet3Pt',
                                           'Jet1Eta', 'Jet2Eta',# 'Jet3Eta',
                                           #'Jet1Phi', 'Jet2Phi', 'Jet3Phi',
                                           'BJet1BDiscr', 'BJet1Pt', 'BJet1Eta', #'BJet1Phi', 
                                           'BJet2BDiscr', 'BJet2Pt', 'BJet2Eta', #'BJet2Phi',
                                           'HT', 'HTs', 'EventBalance', 'Centrality',
                                           'JetMultCharge', 'JetMult', 'BJetMult']

    plotter._variableDict['MET']        = ['Met', 'MHT', 'METLD', 'MHT-MET', 'MetPhi', 'MetSumEt',
                                           'MetLepton1DeltaPhi', 'MetLepton2DeltaPhi'
                                           'MetLepDeltaPhiMin', 'nearLepIndex', 'ProjectedMet', 'MetLepton3DeltaPhi'] 

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

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

        inclusive_plotter.make_save_path(outFile, clean=True)
        p_plot.append(Process(name = cut[2:] + '/inclusive', target = plotter_wrapper, args=(inclusive_plotter, 'inclusive', inFile, outFile, do1D, do2D, doLog, doRatio, doEff)))


    ### 3l selection ###
    if do3l:

        plotter_3l = copy.deepcopy(plotter)
        plotter_3l.add_datasets(samples['3l'], Clear=True)
        plotter_3l._overlayList = ['DATA', 'FCNH']

        for i, cut in enumerate(cutList):
            inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

            if doLog:
                outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
            else:
                outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

            plotter_3l.make_save_path(outFile, clean=True)

            for category in cat3l:
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(plotter_3l, category, inFile, outFile, do1D, do2D, doLog, doRatio, doEff)))

    ### ss selection ###
    if doSS:
        ss_plotter = copy.deepcopy(plotter)
        ss_plotter.add_datasets(samples['ss'], Clear=True)
        ss_plotter._overlayList = ['DATA', 'FCNH']

        for i, cut in enumerate(cutList):
            inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

            if doLog:
                outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
            else:
                outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

            ss_plotter.make_save_path(outFile, clean=True)

            for category in catSS:
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(ss_plotter, category, inFile, outFile, do1D, do2D, doLog, doRatio, doEff)))


    ### os selection ###
    if doOS:
        os_plotter = copy.deepcopy(plotter)
        os_plotter.add_datasets(samples['os'], Clear=True)
        os_plotter._overlayList = ['DATA']

        for i, cut in enumerate(cutList):
            inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

            if doLog:
                outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
            else:
                outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

            os_plotter.make_save_path(outFile, clean=True)

            for category in catOS:
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(os_plotter, category, inFile, outFile, do1D, do2D, doLog, doRatio, doEff)))

    doLog = False

    ### WZ control region
    if 'CR_WZ' in crList:
        wz_plotter = copy.deepcopy(plotter)
        wz_plotter.add_datasets(samples['WZ'], Clear=True)
        wz_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 6, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'CR_WZ')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'CR_WZ')

        wz_plotter.make_save_path(outFile, clean=True)

        for category in cat3l:
            p_plot.append(Process(name = 'CR_WZ/' + category, target = plotter_wrapper, args=(wz_plotter, category, inFile, outFile, do1D, False, doLog, doRatio, doEff)))

    ### ttbar control region
    if 'CR_ttbar' in crList:
        ttbar_plotter = copy.deepcopy(plotter)
        ttbar_plotter.add_datasets(samples['ttbar'],  Clear=True)
        ttbar_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 7, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'CR_ttbar')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'CR_ttbar')

        ttbar_plotter.make_save_path(outFile, clean=True)

        p_plot.append(Process(name = 'CR_ttbar/os_emu', target = plotter_wrapper, args=(ttbar_plotter, 'os_emu', inFile, outFile, do1D, False, doLog, doRatio, doEff)))

    ### ttZ control region
    if 'CR_ttZ' in crList:
        ttZ_plotter = copy.deepcopy(plotter)
        ttZ_plotter.add_datasets(samples['ttZ'],  Clear=True)
        ttZ_plotter._overlayList = ['DATA']

        inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, 8, period, batch)

        if doLog:
            outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, 'CR_ttZ')
        else:
            outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, 'CR_ttZ')

        ttZ_plotter.make_save_path(outFile, clean=True)

        for category in cat3l:
            p_plot.append(Process(name = 'CR_ttZ/' + category, target = plotter_wrapper, args=(ttZ_plotter, category, inFile, outFile, do1D, False, doLog, doRatio, doEff)))


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
        #yieldTable.get_scale_factors()
        yieldTable.get_scale_factors(['FCNH'])

    if do3l:
        #yieldTable._columnList  = samples['3l'] + ['BG', 'DATA', 'FCNH']#, 'Significance'] 
        #yieldTable._columnList  = ['BG', 'DATA', 'FCNC_M125_t', 'FCNC_M125_tbar', 'FCNC_M125_t_semilep', 'FCNC_M125_t_ZZ', 'FCNC_M125_t_TauTau','FCNH']# 'Significance'] 
        yieldTable._columnList  = ['BG', 'DATA', 'FCNH']#, 'Significance'] 

        yieldTable.add_datasets(samples['3l'], Clear = True)
        yieldTable.add_datasets('FCNH')
        yieldTable.add_datasets('DATA')

        yieldTable._rowList = 5*['.'] + ['3 lepton', 'Z removal', 'MET', 'HT', 'b-jet'] + 5*['.'] + ['BDT']

        for category in cat3l:
            yieldTable._category = category
            histDict = yieldTable.get_hist_dict('YieldByCut')
            yieldTable.print_table(histDict, doErrors = True, doEff = False, startBin = 1)

    if doSS:
        #yieldTable._columnList  = samples['ss'] + ['BG', 'DATA', 'FCNH']#, 'Significance'] 
        yieldTable._columnList  = ['BG', 'DATA', 'FCNH']#, 'Significance'] 

        yieldTable.add_datasets(samples['ss'], Clear = True)
        yieldTable.add_datasets('FCNH')
        yieldTable.add_datasets('DATA')

        yieldTable._rowList = ['.', '.', '.', '.', '.','ss lepton', 'Z removal', 'MET', 'HT', 'b-jet']

        for category in catSS:
            yieldTable._category = category
            histDict = yieldTable.get_hist_dict('YieldByCut')
            yieldTable.print_table(histDict, doErrors = True, doEff = False, startBin = 1)

    crCats = {'CR_WZ':'3l_inclusive', 'CR_ttbar':'os_emu', 'CR_ttZ':'3l_inclusive'}
    for i,CR in enumerate(crList):

        yieldTable.set_input_file('fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, i+5, period, batch))
        yieldTable._columnList  = samples[CR[3:]] + ['BG', 'DATA']

        yieldTable.add_datasets(samples[CR[3:]], Clear = True)
        yieldTable.add_datasets('DATA')

        yieldTable._rowList = ['preselection'] + (4+i)*['.'] + [CR[3:]]

        yieldTable._category = crCats[CR]
        histDict = yieldTable.get_hist_dict('YieldByCut')
        yieldTable.print_table(histDict, doErrors = False, doEff = False, startBin = 6)

    ### Special case for ZZ->4l control region ###
    yieldTable.set_input_file('fcncAnalysis/combined_histos/{0}_cut1_{1}_{2}.root'.format(selection, period, batch))
    yieldTable.add_datasets(['ZZ4l', 'DATA'], Clear = True)
    yieldTable._columnList  = ['ZZ4l'] + ['BG', 'DATA']
    yieldTable._rowList = 8*['.'] + ['ZZ4l']

    yieldTable._category = 'inclusive'
    histDict = yieldTable.get_hist_dict('YieldByCut')
    yieldTable.print_table(histDict, doErrors = False, doEff = False, startBin = 6)

    tableFile.close()

    subprocess.call('pdflatex -output-dir=yields yields/yields.tex', shell = True)
    subprocess.call('cp yields/yields.pdf plots/{0}/{1}_{2}_{3}/.'.format(currentDate, selection, batch, suffix), shell = True)
    subprocess.call('cp yields/.yields_tmp.tex plots/{0}/{1}_{2}_{3}/yields.tex'.format(currentDate, selection, batch, suffix), shell = True)
