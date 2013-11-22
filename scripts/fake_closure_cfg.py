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

period      = '2012'
LUMIDATA    = 19.712 

doPlots     = True
doLog       = True
doEff       = False
doRatio     = True
do1D        = True
do2D        = False

doYields    = False

doSS        = True
do3l        = True

### Categories to be plotted ###
catSS       = ['ss_inclusive']
catSS.extend(['ss_mumu', 'ss_ee', 'ss_emu'])
cat3l       = ['3l_inclusive']
cat3l.extend(['3l_OSSF', '3l_SSSF'])
cat3l.extend(['3l_eee', '3l_eemu', '3l_emumu', '3l_mumumu'])

### Samples to be included in stacks ###
samples     = {'all':[], 'inclusive':[], '3l':[], 'ss':[]}

samples['all'].append('higgs')
samples['all'].append('Triboson')
samples['all'].append('ttV')
samples['all'].append('Diboson')
samples['all'].append('top')
samples['all'].append('ZJets')
samples['all'].append('WJets')
samples['all'].append('QCD')
#samples['all'].extend(['ZbbToLL', 'WbbToLNu']) #, 'ZGstar'])

samples['3l'].append('Fakes_Triboson')
samples['3l'].append('Fakes_ttV')
samples['3l'].append('Fakes_ZZ4l')
samples['3l'].append('Fakes_WZJets3LNu')
samples['3l'].append('Fakes_top')
samples['3l'].append('Fakes_ZJets')

#samples['ss'].append('Fakes_Triboson')
#samples['ss'].append('Fakes_ttV')
#samples['ss'].append('Fakes_ZZ4l')
#samples['ss'].append('Fakes_WZJets3LNu')
samples['ss'].append('Fakes_top')
samples['ss'].append('Fakes_WJets')
samples['ss'].append('Fakes_QCD')
samples['ss'].append('Fakes_ZJets')

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
    plotter.get_scale_factors()

    ### VARIABLES ###
    ### First specify the directories in which your
    ### histograms are stored.  If directories are 
    ### not used enter '' as the only entry.  Then 
    ### list all of the variable names you wish to 
    ### plot while giving a key value which is the 
    ### directory that they are located in as a key.

    plotter._directoryList1D            = ['Misc', 'Lepton', 'Lep+Jet', 'Dilepton', 'DileptonOS', 'Trilepton', 'MET', 'Jet']

    plotter._variableDict['Misc']       = ['PvMult', 'YieldByCut', 'YieldByCutRaw', 'EventWeight', 'TriggerStatus', 'BDT'] 

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
                                           'Lepton3BJetDeltaPhi', 'Lepton3BJetDeltaEta', 'Lepton3BJetDeltaR', 'Lepton3BJetDeltaPt',
                                           'Lepton1JetDeltaPhi', 'Lepton1JetDeltaEta', 'Lepton1JetDeltaR', 'Lepton1JetDeltaPt',
                                           'Lepton2JetDeltaPhi', 'Lepton2JetDeltaEta', 'Lepton2JetDeltaR', 'Lepton2JetDeltaPt',
                                           'Lepton3JetDeltaPhi', 'Lepton3JetDeltaEta', 'Lepton3JetDeltaR', 'Lepton3JetDeltaPt',
                                           'DileptonJetDeltaPhi', 'DileptonJetDeltaEta', 'DileptonJetDeltaR', 'DileptonJetDeltaPt',
                                           'DileptonBJetDeltaPhi', 'DileptonBJetDeltaEta', 'DileptonBJetDeltaR', 'DileptonBJetDeltaPt',
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

    ### 3l selection ###
    if do3l:

        plotter_3l = copy.deepcopy(plotter)
        plotter_3l.add_datasets(samples['3l'], Clear=True)
        plotter_3l._overlayList = ['Fakes']

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
        ss_plotter._overlayList = ['Fakes']

        for i, cut in enumerate(cutList):
            inFile  = 'fcncAnalysis/combined_histos/{0}_cut{1}_{2}_{3}.root'.format(selection, str(i+1), period, batch)

            if doLog:
                outFile = 'plots/{0}/{1}_{2}_{3}/log/{4}'.format(currentDate, selection, batch, suffix, cut)
            else:
                outFile = 'plots/{0}/{1}_{2}_{3}/linear/{4}'.format(currentDate, selection, batch, suffix, cut)

            ss_plotter.make_save_path(outFile, clean=True)

            for category in catSS:
                p_plot.append(Process(name = cut[2:] + '/' + category, target = plotter_wrapper, args=(ss_plotter, category, inFile, outFile, do1D, do2D, doLog, doRatio, doEff)))


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